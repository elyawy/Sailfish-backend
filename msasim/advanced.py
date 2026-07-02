"""Advanced simulation models"""

from dataclasses import dataclass

import pathlib
import numpy as np
from typing import List
import warnings

import _Sailfish
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
        self._rng = np.random.default_rng(self.indel_model.get_seed())  # Seed the random number generator for reproducibility
        # initialize simulators for each replacement model
        self.substitution_simulators = []
        for model in self.mixture_model.replacement_models:
            sub_model = SubstitutionModel(model)
            self.substitution_simulators.append(sub_model)

    def simulate(self) -> Msa:
        # Simulation logic for the mixture model
        # Step 1: Simulate the indel model to get the template MSA
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            indel_simulator = Simulator(simProtocol=self.indel_model)
            
        sim_context = indel_simulator.protocol.get_sim_context()
        template_msa = indel_simulator()
        sparse_msa: List[List[int]] = template_msa.get_sparse_msa()
        

        msa_length  = template_msa.get_length()

        lengths_per_model =self._rng.multinomial(msa_length, self.mixture_model.model_weights)
        model_ids = np.repeat(np.arange(len(self.mixture_model.replacement_models)), lengths_per_model)
        perm = self._rng.permutation(msa_length)
        model_ids = model_ids[perm]

        substitution_only_msas = []
        
        for model_idx, length in enumerate(lengths_per_model):
            # Step 2: Sample a replacement model based on the weights
            current_model: SubstitutionModel = self.substitution_simulators[model_idx]
            # Step 3: Simulate a single site using the sampled replacement model
            sub_simulator = current_model.build_substitution_simulator(sim_context)
            temp_msa = Msa(length, sim_context)

            sequences = sub_simulator.simulate_substitutions(length, "" ,temp_msa.get_root_positions_in_msa())

            substitution_only_msas.append(sequences)
        
        # concat the simulated sequences from each model to form the final MSA
        # after a row concat shuffle with the model_ids to ensure the sites are in the correct order
        # underlying object in sequences is SparseSequenceContainer which is a vector of SparseSequence objects,
        # each of which is a std::string on the C++ side.
        # sparse_msa is a std::vector<std::vector<int>> encoding.
        # 


        num_sequences = template_msa.get_num_sequences()
        sparse_rows = []

        for row in range(num_sequences):
            # Step 6: concat this row's chars across models, grouped order
            grouped_seq = "".join(substitution_only_msas[m][row] for m in range(len(lengths_per_model)))

            # Step 7: reorder into final site order using perm
            grouped_arr = np.array(list(grouped_seq))
            final_arr = grouped_arr[perm]

            # Step 8: strip gaps using this row's run-length blocks
            sparse_chars = []
            pos = 0
            for block in sparse_msa[row]:
                if block < 0:
                    pos += -block  # gap block, skip in dense array
                else:
                    sparse_chars.append("".join(final_arr[pos:pos + block]))
                    pos += block
            sparse_rows.append("".join(sparse_chars))

        # Step 9
        template_msa.fill_substitutions(_Sailfish.SparseSequenceContainer(sparse_rows))

        return template_msa