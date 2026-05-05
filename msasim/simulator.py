"""MSA simulator"""

import _Sailfish
import warnings
import pathlib
from typing import Optional, List
from .protocol import SimProtocol
from .msa import Msa
from .constants import SIMULATION_TYPE, SIMULATION_TYPES, DNA_MODELS, PROTEIN_MODELS
from .substitutions import SubstitutionModel


class Simulator:
    """Simulate MSAs based on SimProtocol"""

    def __init__(
        self,
        simProtocol: Optional[SimProtocol] = None,
        simulation_type: Optional[SIMULATION_TYPE] = None,
    ):
        if simProtocol is None:
            simProtocol = SimProtocol.default()

        if simProtocol._verify_sim_protocol():
            self._simProtocol = simProtocol
            self._indel_simulator = _Sailfish.IndelSimulator(
                self._simProtocol.get_sim_context(), self._simProtocol._sim_protocol
            )
        else:
            raise ValueError("failed to verify simProtocol")

        if not simulation_type:
            warnings.warn("simulation type not provided -> running indel only simulation")
            simulation_type = SIMULATION_TYPE.NOSUBS
            self._simProtocol._sim_protocol.set_site_rate_model(_Sailfish.SiteRateModel.SIMPLE)

        if simulation_type not in SIMULATION_TYPES:
            raise ValueError(
                f"unknown simulation type, please provide one of the following: "
                f"{[e.name for e in SIMULATION_TYPE]}"
            )

        self._simulation_type = simulation_type
        if self._simulation_type == SIMULATION_TYPE.NOSUBS:
            self._substitution_simulator = None
            self._sub_model = None
        else:
            self._sub_model = SubstitutionModel(simulation_type)
            sim_context = self._simProtocol.get_sim_context()
            if self._simulation_type == SIMULATION_TYPE.PROTEIN:
                self._substitution_simulator = _Sailfish.AminoSubstitutionSimulator(
                    self._sub_model.factory, sim_context
                )
            else:
                self._substitution_simulator = _Sailfish.NucleotideSubstitutionSimulator(
                    self._sub_model.factory, sim_context
                )

        # Compute flags once — used to pick strategy at __init__ time
        self._has_indels = not (
            self._simProtocol._is_insertion_rate_zero
            and self._simProtocol._is_deletion_rate_zero
        )
        self._has_subs = self._simulation_type != SIMULATION_TYPE.NOSUBS
        self._is_indel_aware = (
            self._simProtocol.get_site_rate_model() == _Sailfish.SiteRateModel.INDEL_AWARE
        )

        # Pick strategy once — no branching inside simulate()
        if not self._has_indels and not self._has_subs:
            self._strategy = self._simulate_root_only
        elif not self._has_subs:
            self._strategy = self._simulate_indels_only
        elif not self._has_indels:
            self._strategy = self._simulate_subs_only
        else:
            self._strategy = self._simulate_full

    # ------------------------------------------------------------------
    # Private simulation strategies
    # ------------------------------------------------------------------

    def _build_msa_no_indels(self) -> Msa:
        """Build a trivial MSA directly from root sequence size (no indel events)."""
        return Msa(self._simProtocol.get_sequence_size(), self._simProtocol.get_sim_context())

    def _build_msa_with_indels(self) -> Msa:
        """Run indel simulation and build MSA, setting up the category sampler if INDEL_AWARE."""
        sim_context = self._simProtocol.get_sim_context()
        eventmap = self.generate_events()
        if self._is_indel_aware:
            # INDEL_AWARE: the MSA builder assigns rate categories to inserted sites
            # during construction, so the sampler must be ready before Msa() is called.
            category_sampler = self._sub_model.factory.get_rate_category_sampler(
                self._simProtocol.get_max_insertion_length()
            )
            sim_context.set_category_sampler(category_sampler)
        return Msa(eventmap, sim_context)

    def _apply_substitutions(self, msa: Msa) -> None:
        """Run substitution simulation and fill the MSA in-place."""
        rate_categories = msa.get_per_site_rate_categories()
        self._substitution_simulator.set_per_site_rate_categories(rate_categories)
        self._substitution_simulator.set_aligned_sequence_map(msa._msa)
        substitutions = self._substitution_simulator.simulate_substitutions(msa.get_length())
        msa.fill_substitutions(substitutions)

    def _apply_substitutions_to_disk(self, msa: Msa, output_path: pathlib.Path) -> None:
        """Simulate substitutions and write directly to disk without holding full MSA in memory."""
        self._substitution_simulator.set_aligned_sequence_map(msa._msa)
        self._substitution_simulator.set_per_site_rate_categories(
            msa.get_per_site_rate_categories()
        )
        self._substitution_simulator.simulate_and_write_substitutions(
            msa.get_length(), str(output_path)
        )

    def _simulate_root_only(self, output_path: Optional[pathlib.Path]) -> Optional[Msa]:
        """No indels, no substitutions — return the bare root sequence MSA."""
        msa = self._build_msa_no_indels()
        if output_path:
            msa.write_msa(str(output_path))
            return None
        return msa

    def _simulate_indels_only(self, output_path: Optional[pathlib.Path]) -> Optional[Msa]:
        """Indels only, no substitutions — template alignment with gap structure."""
        msa = self._build_msa_with_indels()
        if output_path:
            msa.write_msa(str(output_path))
            return None
        return msa

    def _simulate_subs_only(self, output_path: Optional[pathlib.Path]) -> Optional[Msa]:
        """Substitutions only, no indels — linear MSA with no gap structure."""
        msa = self._build_msa_no_indels()
        if output_path:
            self._apply_substitutions_to_disk(msa, output_path)
            return None
        self._apply_substitutions(msa)
        return msa

    def _simulate_full(self, output_path: Optional[pathlib.Path]) -> Optional[Msa]:
        """Full simulation: indels + substitutions."""
        msa = self._build_msa_with_indels()
        if output_path:
            self._apply_substitutions_to_disk(msa, output_path)
            return None
        self._apply_substitutions(msa)
        return msa

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def simulate(
        self, output_path: Optional[pathlib.Path] = None, times: int = 1, 
    ) -> Optional[List[Msa]]:
        """
        Run the simulation.

        Args:
            times:       Number of replicates to generate. Ignored when output_path is set.
            output_path: If provided, write directly to disk (low-memory mode) and return None.
                         Otherwise collect and return a list of Msa objects.
        """
        if output_path is not None:
            self._strategy(output_path)
            return [None]  # Return a list of None for consistency with return type
        

        return [self._strategy(None) for _ in range(times)]

    def __call__(self, output_path: Optional[pathlib.Path] = None) -> Msa:
        return self.simulate(output_path=output_path, times=1)[0]

    # ------------------------------------------------------------------
    # Configuration / accessors (unchanged)
    # ------------------------------------------------------------------

    def reset_substitution_simulator(self, modelFactory: _Sailfish.modelFactory) -> None:
        if self._simulation_type == SIMULATION_TYPE.PROTEIN:
            self._substitution_simulator = _Sailfish.AminoSubstitutionSimulator(
                modelFactory, self._simProtocol.get_sim_context()
            )
        else:
            self._substitution_simulator = _Sailfish.NucleotideSubstitutionSimulator(
                modelFactory, self._simProtocol.get_sim_context()
            )

    def set_replacement_model(
        self,
        model: _Sailfish.modelCode,
        amino_model_file: pathlib.Path = None,
        model_parameters: List = None,
        gamma_parameters_alpha: float = 1.0,
        gamma_parameters_categories: int = 1,
        invariant_sites_proportion: float = 0.0,
        site_rate_correlation: float = 0.0,
    ) -> None:
        next_simulation_type = (
            SIMULATION_TYPE.PROTEIN if model in PROTEIN_MODELS else SIMULATION_TYPE.DNA
        )
        self._sub_model.set_replacement_model(
            model=model,
            amino_model_file=amino_model_file,
            model_parameters=model_parameters,
            gamma_parameters_alpha=gamma_parameters_alpha,
            gamma_parameters_categories=gamma_parameters_categories,
            invariant_sites_proportion=invariant_sites_proportion,
            site_rate_correlation=site_rate_correlation,
            simulation_type=self._simulation_type,
        )
        if next_simulation_type != self._simulation_type:
            raise ValueError(
                f"replacement model {model} is not compatible with current simulation type "
                f"{self._simulation_type.name}. Please initialize a separate Simulator."
            )
        else:
            self._substitution_simulator.init_substitution_sim(self._sub_model.factory)

    @property
    def protocol(self) -> SimProtocol:
        return self._simProtocol

    def get_sequences_to_save(self) -> List[bool]:
        return self._simProtocol.get_sim_context().get_nodes_to_save()

    def save_root_sequence(self):
        self._simProtocol.get_sim_context().set_save_root()

    def save_all_nodes_sequences(self):
        self._simProtocol.get_sim_context().set_save_all_nodes()

    def save_leaves_sequences(self):
        self._simProtocol.get_sim_context().set_save_leaves()

    def generate_events(self) -> List[List[_Sailfish.IndelEvent]]:
        return self._indel_simulator.generate_events()

    def save_rates(self, is_save: bool) -> None:
        self._substitution_simulator.set_save_rates(is_save)

    def get_rates(self) -> List[float]:
        return self._substitution_simulator.get_site_rates()

    def get_rate_categories(self) -> List[int]:
        return self._substitution_simulator.get_per_site_rate_categories()