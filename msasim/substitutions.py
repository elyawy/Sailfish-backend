"""Substitution model configuration and rate model management."""

import _Sailfish
import warnings
import pathlib
from typing import List, Optional

from .constants import MODEL_CODES, SIMULATION_TYPE


# Default model applied at construction — cheap to build, valid for both DNA and protein-less sims
_DEFAULT_MODEL = {SIMULATION_TYPE.DNA: MODEL_CODES.NUCJC, SIMULATION_TYPE.PROTEIN: MODEL_CODES.JONES}
_DEFAULT_GAMMA_ALPHA = 1.0
_DEFAULT_GAMMA_CATEGORIES = 1


class SubstitutionModel:
    """
    Owns the C++ modelFactory and all substitution/rate model configuration.

    Tracks the substitution model identity (model code, parameters, amino file)
    so that the expensive replacement model rebuild is skipped when only
    site-rate parameters (gamma, invariants, correlation) change.
    """

    def __init__(self, model_type: SIMULATION_TYPE) -> None:
        self._factory = _Sailfish.modelFactory()

        # Track substitution model identity to detect real changes
        self._current_model: Optional[MODEL_CODES] = None
        self._current_model_parameters: Optional[List] = None
        self._current_amino_file: Optional[str] = None
        self._model_type: SIMULATION_TYPE = model_type
        # Apply defaults so the factory is immediately in a valid COMPLETE state
        self._apply_default_model()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    @property
    def factory(self) -> _Sailfish.modelFactory:
        """The underlying C++ modelFactory — passed to SubstitutionSimulator."""
        return self._factory

    def set_replacement_model(
        self,
        model: _Sailfish.modelCode,
        amino_model_file: pathlib.Path = None,
        model_parameters: List = None,
        gamma_parameters_alpha: float = 1.0,
        gamma_parameters_categories: int = 1,
        invariant_sites_proportion: float = 0.0,
        site_rate_correlation: float = 0.0,
        simulation_type: SIMULATION_TYPE = None,
    ) -> None:
        """
        Configure the substitution and site-rate model.

        The C++ replacement model (expensive) is only rebuilt when the substitution
        model itself changes. Updating only site-rate parameters (gamma, invariants,
        correlation) reuses the cached model.
        """
        if not simulation_type:
            raise ValueError("simulation_type is required to set the replacement model")

        if not model:
            raise ValueError(
                f"please provide a substitution model from the following list: {_Sailfish.modelCode}"
            )
        if int(gamma_parameters_categories) != gamma_parameters_categories:
            raise ValueError(
                f"gamma_parameters_categories has to be a positive int value, "
                f"received: {gamma_parameters_categories}"
            )

        # Validate per simulation type
        if simulation_type == SIMULATION_TYPE.PROTEIN:
            if model_parameters:
                raise ValueError(
                    f"no model parameters are used in protein models, "
                    f"received: {model_parameters}"
                )
        else:
            if model == MODEL_CODES.NUCJC and model_parameters:
                raise ValueError("no model parameters in JC model, received: {model_parameters}")
            if model != MODEL_CODES.NUCJC and not model_parameters:
                raise ValueError("please provide model_parameters for this nucleotide model")

        sub_model_changed = (
            model != self._current_model
            or model_parameters != self._current_model_parameters
            or (str(amino_model_file) if amino_model_file else None) != self._current_amino_file
        )

        if sub_model_changed:
            self._factory.reset()
            self._factory.set_replacement_model(model)

            if simulation_type == SIMULATION_TYPE.PROTEIN:
                if model == MODEL_CODES.CUSTOM and amino_model_file:
                    self._factory.set_amino_replacement_model_file(str(amino_model_file))
            else:
                if model_parameters:
                    self._factory.set_model_parameters(model_parameters)

            self._current_model = model
            self._current_model_parameters = model_parameters
            self._current_amino_file = str(amino_model_file) if amino_model_file else None

        # Always update site-rate model — cheap, does not touch the cached pij
        rates, probs, transition_matrix = self._create_site_rate_model(
            gamma_alpha=gamma_parameters_alpha,
            gamma_categories=gamma_parameters_categories,
            invariant_proportion=invariant_sites_proportion,
            site_rate_correlation=site_rate_correlation,
        )
        self._factory.set_site_rate_model(rates, probs, transition_matrix)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _apply_default_model(self) -> None:
        """Set factory to a valid COMPLETE state using NUCJC + single rate category."""
        self._factory.set_replacement_model(_DEFAULT_MODEL[self._model_type])
        rates, probs, transition_matrix = self._create_site_rate_model(
            gamma_alpha=_DEFAULT_GAMMA_ALPHA,
            gamma_categories=_DEFAULT_GAMMA_CATEGORIES,
        )
        self._factory.set_site_rate_model(rates, probs, transition_matrix)
        self._current_model = _DEFAULT_MODEL[self._model_type]
        self._current_model_parameters = None
        self._current_amino_file = None

    def _create_site_rate_model(
        self,
        gamma_alpha: float = 1.0,
        gamma_categories: int = 1,
        invariant_proportion: float = 0.0,
        site_rate_correlation: float = 0.0,
    ) -> tuple:
        """
        Compute rate categories, stationary probabilities, and transition matrix.

        Returns:
            (rates, probs, transition_matrix)
        """
        if invariant_proportion < 0.0 or invariant_proportion >= 1.0:
            raise ValueError(
                f"invariant_proportion must be in [0, 1), received: {invariant_proportion}"
            )
        if site_rate_correlation < 0.0 or site_rate_correlation >= 1.0:
            raise ValueError(
                f"site_rate_correlation must be in [0, 1), received: {site_rate_correlation}"
            )
        if site_rate_correlation > 0.0 and gamma_categories == 1:
            warnings.warn(
                "site_rate_correlation > 0 requires gamma_categories > 1. "
                "Setting site_rate_correlation to 0.0"
            )
            site_rate_correlation = 0.0

        gamma_dist = _Sailfish.GammaDistribution(gamma_alpha, gamma_categories)
        rates = list(gamma_dist.getAllRates())
        probs = list(gamma_dist.getAllRatesProb())

        if invariant_proportion > 0.0:
            scale_factor = 1.0 - invariant_proportion
            probs = [p * scale_factor for p in probs]
            rates.insert(0, 0.0)
            probs.insert(0, invariant_proportion)

        transition_matrix = []
        if site_rate_correlation > 0.0:
            if invariant_proportion > 0.0:
                warnings.warn(
                    "site_rate_correlation and invariant_sites_proportion cannot be used together. "
                    "Using invariant sites only, ignoring correlation."
                )
            else:
                try:
                    from msasim.correlation import build_auto_gamma_transition_matrix
                    transition_matrix = build_auto_gamma_transition_matrix(
                        alpha=gamma_alpha,
                        categories=gamma_categories,
                        rho=site_rate_correlation,
                    )
                except ImportError:
                    warnings.warn(
                        "site_rate_correlation > 0 requires scipy. "
                        "Install with: pip install scipy or pip install 'msasim[correlation]'. "
                        "Ignoring correlation parameter."
                    )

        if len(transition_matrix) == 0:
            transition_matrix = [probs for _ in range(len(probs))]

        return rates, probs, transition_matrix
