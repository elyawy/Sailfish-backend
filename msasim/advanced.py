"""Advanced simulation models"""

from dataclasses import dataclass

import pathlib
import numpy as np
from typing import List
from .protocol import SimProtocol
from .msa import Msa
from .substitutions import ReplacementModelSpec, SubstitutionModel
from .simulator import Simulator

@dataclass
class PartitionModel:
    """Partition model for advanced simulations.

    Attributes:
        name (str): Name of the partition model.
        replacement_model (ReplacementModelSpec): Substitution model used in the partition.
        indel_model (SimProtocol): Indel model used in the partition.
        """
    name: str
    replacement_model: ReplacementModelSpec
    indel_model: SimProtocol



class Partitions:
    def __init__(self):
        self._partitions: List[PartitionModel] = []
        

    def add_partition(self, partition: PartitionModel):
        self._partitions.append(partition)

    def merge_msas_to_file(self, msas: List[Msa], output_path: pathlib.Path) -> None:
        partition_dicts = []
        for msa in msas:
            rows = {}
            for i in range(msa.get_num_sequences()):
                name, seq = msa.get_msa_row(i).split("\n", 1)
                rows[name.lstrip(">")] = seq
            partition_dicts.append((rows, msa.get_length()))

        all_taxa = set()
        for rows, _ in partition_dicts:
            all_taxa.update(rows.keys())

        with open(output_path, "w") as f:
            for taxon in sorted(all_taxa):
                sequence = "".join(
                    rows.get(taxon, "-" * length) for rows, length in partition_dicts
                )
                f.write(f">{taxon}\n{sequence}\n")
                
        
            
    
    def simulate(self, output_path: pathlib.Path) -> None:
        if not self._partitions:
            raise ValueError("no partitions added, use add_partition() first")
        msa_results = []
        for partition in self._partitions:
            # Simulation logic for each partition
            simulator = Simulator(simProtocol=partition.indel_model,
                                  replacement_model=partition.replacement_model)
            
            msa = simulator()
            msa_results.append(msa)

        self.merge_msas_to_file(msa_results, output_path)


@dataclass
class MixtureModel:
    """Mixture model for advanced simulations.

    Attributes:
        name (str): Name of the mixture model.
        replacement_models (List[ReplacementModelSpec]): Substitution models used in the mixture.
        model_weights (List[float]): Weights for each substitution model, must sum to 1.0.
        indel_model (SimProtocol): Indel model used in the mixture.
        """
    name: str
    replacement_models: List[ReplacementModelSpec]
    model_weights: List[float] # must sum to 1.0
    indel_model: SimProtocol

    def __post_init__(self):
        # all models must be of the same type (nucleotide or amino acid)
        model_types = {model.model_type for model in self.replacement_models}
        if len(model_types) > 1:
            raise ValueError("All replacement models must be of the same type (nucleotide or amino acid)")
        if len(self.replacement_models) != len(self.model_weights):
            raise ValueError("replacement_models and model_weights must have the same length")
        if not (0.999 <= sum(self.model_weights) <= 1.001):  # Allowing for floating point error
            raise ValueError("model_weights must sum to 1.0")

class Mixture:
    """
    Mixture simulation model for advanced simulations.
    How this works: We intialize separately a single indle only simulator, 
    and a substitution only simulator for each of the models in the mixture.
    the indel only sim creates the template msa and decides the length of the overall msa.
    the substitution only sims are then sampled at random for each site - each one then simulates a single site
    the sites are combined and the MSA is filled with them.
    finally we return the MSA to the user.
    """
    def __init__(self, mixture_model: MixtureModel):
        self.mixture_model = mixture_model
        self.indel_model = mixture_model.indel_model
        # initialize simulators for each replacement model
        self.substitution_simulators = []
        for model in self.mixture_model.replacement_models:
            sub_model = SubstitutionModel(model)
            self.substitution_simulators.append(sub_model)

    def simulate(self) -> None:
        # Simulation logic for the mixture model
        # Step 1: Simulate the indel model to get the template MSA
        indel_simulator = Simulator(simProtocol=self.indel_model)
        template_msa = indel_simulator()
        sparse_msa = template_msa.get_sparse_msa()
        

        msa_length  = template_msa.get_length()

        np.random.multinomial()
        
        for i in range(msa_length):
            # Step 2: Sample a replacement model based on the weights
            sampled_model: SubstitutionModel = random.choices(
                self.substitution_simulators, 
                weights=self.mixture_model.model_weights, 
                k=1
            )[0]
            # Step 3: Simulate a single site using the sampled replacement model
            sub_simulator = sampled_model.build_substitution_simulator()
            sequences = sub_simulator.simulate_substitutions(1)  # Simulate a single site

            # Step 4: Fill the template MSA with the simulated site
            template_msa.fill_site(i, site_msa)
