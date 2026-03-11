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
            simProtocol = SimProtocol.default()  # fresh object each time

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

            # Initialize the C++ substitution simulator immediately with the default model
            sim_context = self._simProtocol.get_sim_context()
            if self._simulation_type == SIMULATION_TYPE.PROTEIN:
                self._substitution_simulator = _Sailfish.AminoSubstitutionSimulator(
                    self._sub_model.factory, sim_context
                )
            else:
                self._substitution_simulator = _Sailfish.NucleotideSubstitutionSimulator(
                    self._sub_model.factory, sim_context
                )

    def reset_substitution_simulator(self, modelFactory: _Sailfish.modelFactory) -> None:
        if self._simulation_type == SIMULATION_TYPE.PROTEIN:
            self._substitution_simulator = _Sailfish.AminoSubstitutionSimulator(modelFactory, self._simProtocol.get_sim_context())
        else:
            self._substitution_simulator = _Sailfish.NucleotideSubstitutionSimulator(modelFactory, self._simProtocol.get_sim_context())

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
        next_simulation_type = SIMULATION_TYPE.PROTEIN if model in PROTEIN_MODELS else SIMULATION_TYPE.DNA
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
                f"replacement model {model} is not compatible with current simulation type {self._simulation_type.name}. Please initialize a separate Simulator with the desired simulation type."
            )
        else:
            self._substitution_simulator.init_substitution_sim(self._sub_model.factory)

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

    def gen_substitutions(self, msa: Msa):
        rate_categories = msa.get_per_site_rate_categories()
        self._substitution_simulator.set_per_site_rate_categories(rate_categories)
        self._substitution_simulator.set_aligned_sequence_map(msa._msa)
        return self._substitution_simulator.simulate_substitutions(msa.get_length())

    def simulate(self, times: int = 1) -> List[Msa]:
        Msas = []
        sim_context = self._simProtocol.get_sim_context()
        for _ in range(times):
            if self._simProtocol._is_insertion_rate_zero and self._simProtocol._is_deletion_rate_zero:
                msa = Msa(self._simProtocol.get_sequence_size(), sim_context)
            else:                    
                eventmap = self.generate_events()
                if self._simulation_type != SIMULATION_TYPE.NOSUBS:
                    category_sampler = self._sub_model.factory.get_rate_category_sampler(
                        self._simProtocol.get_max_insertion_length()
                    )
                    sim_context.set_category_sampler(category_sampler)
                msa = Msa(eventmap, sim_context)

            if self._simulation_type != SIMULATION_TYPE.NOSUBS:
                substitutions = self.gen_substitutions(msa)
                msa.fill_substitutions(substitutions)

            Msas.append(msa)
        return Msas

    def simulate_low_memory(self, output_file_path: pathlib.Path) -> None:
        sim_context = self._simProtocol.get_sim_context()
        if self._simProtocol._is_insertion_rate_zero and self._simProtocol._is_deletion_rate_zero:
            msa_length = self._simProtocol.get_sequence_size()
            msa = Msa(msa_length, sim_context)
        else:
            eventmap = self.generate_events()
            category_sampler = self._sub_model.factory.get_rate_category_sampler(
                self._simProtocol.get_max_insertion_length()
            )
            sim_context.set_category_sampler(category_sampler)
            msa = Msa(eventmap, sim_context)
            msa_length = msa.get_length()
            self._substitution_simulator.set_aligned_sequence_map(msa._msa)
            self._substitution_simulator.set_per_site_rate_categories(
                msa.get_per_site_rate_categories()
            )

        if self._simulation_type == SIMULATION_TYPE.NOSUBS:
            msa.write_msa(str(output_file_path))
        else:
            self._substitution_simulator.simulate_and_write_substitutions(
                msa_length, str(output_file_path)
            )

    def __call__(self) -> Msa:
        return self.simulate(1)[0]

    def save_rates(self, is_save: bool) -> None:
        self._substitution_simulator.set_save_rates(is_save)

    def get_rates(self) -> List[float]:
        return self._substitution_simulator.get_site_rates()

    def get_rate_categories(self) -> List[int]:
        return self._substitution_simulator.get_per_site_rate_categories()