"""Substitution model configuration and rate model management."""

from dataclasses import dataclass

import _Sailfish
import warnings
import pathlib
from typing import List, Optional

from .constants import MODEL_CODES, ALPHABET_CODES, ALPHABET_SIZES


_DEFAULT_GAMMA_ALPHA = 1.0
_DEFAULT_GAMMA_CATEGORIES = 1

@dataclass
class SiteRateModelSpec:
    """Site rate model specification for substitution models.

    Supports two modes (mutually exclusive):
    - Gamma: discretized gamma distribution with optional invariant sites and site rate correlation.
    - Free rates: explicit rates and weights, summing to 1.0.

    Attributes:
        gamma_alpha (float): Shape parameter for the gamma distribution.
        gamma_categories (int): Number of discrete rate categories.
        invariant_proportion (float): Proportion of invariant sites.
        site_rate_correlation (float): Autocorrelation between adjacent site rates.
        free_rates (Optional[List[float]]): Explicit rate values (free rates mode).
        free_rate_weights (Optional[List[float]]): Weights for each free rate, must sum to 1.0.
    """
    gamma_alpha: float = _DEFAULT_GAMMA_ALPHA
    gamma_categories: int = _DEFAULT_GAMMA_CATEGORIES
    invariant_proportion: float = 0.0
    site_rate_correlation: float = 0.0
    free_rates: Optional[List[float]] = None
    free_rate_weights: Optional[List[float]] = None

    # Validate on creation
    def __post_init__(self):
        # free rates take precedence over gamma model if both are provided
        if (self.free_rates is None) != (self.free_rate_weights is None):
            raise ValueError("free_rates and free_rate_weights must both be set or both be None")

        if self.free_rates is not None and self.free_rate_weights is not None:
            if len(self.free_rates) != len(self.free_rate_weights):
                raise ValueError(
                    f"free_rates and free_rate_weights must have the same length, "
                    f"received: len(free_rates)={len(self.free_rates)}, len(free_rate_weights)={len(self.free_rate_weights)}"
                )
            if any(rate <= 0.0 for rate in self.free_rates):
                raise ValueError(
                    f"all free_rates must be positive, received: {self.free_rates}"
                )
            if any(weight < 0.0 for weight in self.free_rate_weights):
                raise ValueError(
                    f"all free_rate_weights must be non-negative, received: {self.free_rate_weights}"
                )
            total_weight = sum(self.free_rate_weights)
            if abs(total_weight - 1.0) > 1e-6:
                raise ValueError(
                    f"sum of free_rate_weights must be 1.0, received: {self.free_rate_weights}"
                )
        
        if self.gamma_categories <= 0:
            raise ValueError(
                f"gamma_categories must be a positive integer, received: {self.gamma_categories}"
            )
        if self.gamma_alpha <= 0.0:
            raise ValueError(
                f"gamma_alpha must be a positive float, received: {self.gamma_alpha}"
            )
        if self.invariant_proportion < 0.0 or self.invariant_proportion >= 1.0:
            raise ValueError(
                f"invariant_proportion must be in [0, 1), received: {self.invariant_proportion}"
            )
        if self.site_rate_correlation < 0.0 or self.site_rate_correlation >= 1.0:
            raise ValueError(
                f"site_rate_correlation must be in [0, 1), received: {self.site_rate_correlation}"
            )
        if self.site_rate_correlation > 0.0 and self.gamma_categories == 1:
            raise ValueError(
                "site_rate_correlation > 0 requires gamma_categories > 1. "
                f"received: site_rate_correlation={self.site_rate_correlation}, gamma_categories={self.gamma_categories}"
            )
    

@dataclass
class ReplacementModelSpec:
    """Replacement model specification for substitution models.

    Attributes:
        model (MODEL_CODES): Substitution model code.
        alphabet (ALPHABET_CODES): Alphabet type (DNA, protein, codon).
        amino_model_file (Optional[pathlib.Path]): Path to amino acid model file (for protein models).
        model_parameters (Optional[List]): List of model parameters (for nucleotide models).
    """
    model: MODEL_CODES
    alphabet: ALPHABET_CODES
    model_parameters: Optional[List] = None
    amino_model_file: Optional[pathlib.Path] = None
    site_rate_model: Optional[SiteRateModelSpec] = None

        
    #validate on creation
    def __post_init__(self):

        # Validate per simulation type
        if self.alphabet == ALPHABET_CODES.PROTEIN:
            if self.model_parameters and (self.model != MODEL_CODES.NONREVERSIBLE):
                raise ValueError(
                    f"no model parameters are used in protein models, "
                    f"received: {self.model_parameters}"
                )
        elif self.alphabet == ALPHABET_CODES.CODON:
            if self.model == MODEL_CODES.EMPIRICODON and self.model_parameters:
                raise ValueError(
                    f"no model parameters in EMPIRICODON model, received: {self.model_parameters}"
                )
            if self.model == MODEL_CODES.WYANG and not self.model_parameters:
                raise ValueError(
                    f"Selection factor and (Transition/Transversion) ratio parameters are required in WYANG model, received: {self.model_parameters}"
                )
        else:
            if self.model == MODEL_CODES.NUCJC and self.model_parameters:
                raise ValueError("no model parameters in JC model, received: {self.model_parameters}")
            if self.model != MODEL_CODES.NUCJC and not self.model_parameters:
                raise ValueError("please provide model_parameters for this nucleotide model")

        if self.model == MODEL_CODES.REVERSIBLE and not self.amino_model_file:
            raise ValueError("amino_model_file is required for custom REVERSIBLE protein model")
        
        if self.model == MODEL_CODES.NONREVERSIBLE:
            alphabet_size = ALPHABET_SIZES[self.alphabet]
            if self.model_parameters is None or len(self.model_parameters) != (alphabet_size * alphabet_size + alphabet_size):
                raise ValueError("model_parameters filled with Q matrix and frequencies are required for NONREVERSIBLE models")

class SubstitutionModel:
    """
    Owns the C++ modelFactory and all substitution/rate model configuration.

    Tracks the substitution model identity (model code, parameters, amino file)
    so that the expensive replacement model rebuild is skipped when only
    site-rate parameters (gamma, invariants, correlation) change.
    """

    def __init__(self, replacement_model_spec: ReplacementModelSpec) -> None:
        self._factory = _Sailfish.modelFactory()
        self._spec: Optional[ReplacementModelSpec] = None
        self.set_replacement_model(replacement_model_spec)
    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    @property
    def factory(self) -> _Sailfish.modelFactory:
        """The underlying C++ modelFactory — passed to SubstitutionSimulator."""
        return self._factory
    
    @property
    def alphabet(self) -> ALPHABET_CODES:
        """The simulation type (DNA or protein) based on the current replacement model."""
        return self._spec.alphabet

    def set_replacement_model(
        self,
        replacement_model_spec: ReplacementModelSpec,
    ) -> None:
        """
        Configure the substitution and site-rate model.

        The C++ replacement model (expensive) is only rebuilt when the substitution
        model itself changes. Updating only site-rate parameters (gamma, invariants,
        correlation) reuses the cached model.
        """

        sub_model_changed = (
            self._spec is None
            or replacement_model_spec.model != self._spec.model
            or replacement_model_spec.model_parameters != self._spec.model_parameters
            or replacement_model_spec.amino_model_file != self._spec.amino_model_file
        )


        if sub_model_changed:
            self._factory.build_model(replacement_model_spec.alphabet, 
                                      replacement_model_spec.model, 
                                      replacement_model_spec.model_parameters or [],
                                      str(replacement_model_spec.amino_model_file) if replacement_model_spec.amino_model_file else "",)
 
        # Always update site-rate model — cheap, does not touch the cached pij
        rates, probs, transition_matrix = self._create_site_rate_model(
            replacement_model_spec.site_rate_model or SiteRateModelSpec()
        )
        self._factory.set_site_rate_model(rates, probs, transition_matrix)
 
        self._spec = replacement_model_spec
    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _create_site_rate_model(self, site_rate_model: SiteRateModelSpec) -> tuple:
        """
        Compute rate categories, stationary probabilities, and transition matrix.
 
        Returns:
            (rates, probs, transition_matrix)
        """
        if site_rate_model.free_rates is not None:
            rates = list(site_rate_model.free_rates)
            probs = list(site_rate_model.free_rate_weights)
        else:
            gamma_dist = _Sailfish.GammaDistribution(
                site_rate_model.gamma_alpha, site_rate_model.gamma_categories
            )
            rates = list(gamma_dist.getAllRates())
            probs = list(gamma_dist.getAllRatesProb())
 
            if site_rate_model.invariant_proportion > 0.0:
                scale_factor = 1.0 - site_rate_model.invariant_proportion
                probs = [p * scale_factor for p in probs]
                rates.insert(0, 0.0)
                probs.insert(0, site_rate_model.invariant_proportion)
 
        transition_matrix = []
        if site_rate_model.site_rate_correlation > 0.0:
            if site_rate_model.invariant_proportion > 0.0:
                warnings.warn(
                    "site_rate_correlation and invariant_proportion cannot be used together. "
                    "Using invariant sites only, ignoring correlation."
                )
            else:
                try:
                    from msasim.correlation import build_auto_gamma_transition_matrix
                    transition_matrix = build_auto_gamma_transition_matrix(
                        categories=site_rate_model.gamma_categories,
                        rho=site_rate_model.site_rate_correlation,
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
 
    def build_substitution_simulator(self, sim_context):
        if self.model_type == ALPHABET_CODES.PROTEIN:
            cls = _Sailfish.AminoSubstitutionSimulator
        elif self.model_type == ALPHABET_CODES.DNA:
            cls = _Sailfish.NucleotideSubstitutionSimulator
        elif self.model_type == ALPHABET_CODES.CODON:
            cls = _Sailfish.CodonSubstitutionSimulator
        return cls(self.factory, sim_context)
