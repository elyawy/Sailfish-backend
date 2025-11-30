"""MSA simulator"""

import _Sailfish
import warnings
import pathlib
from typing import Dict, Optional, List
from .protocol import SimProtocol
from .distributions import PoissonDistribution
from .msa import Msa
from .constants import MODEL_CODES, SIMULATION_TYPE


# TODO delete one of this (I think the above if not used)
class BlockTreePython:
    '''
    Used to contain the events on a multiple branches (entire tree).
    '''
    def __init__(self, branch_block_dict: Dict[str, _Sailfish.Block]):
        self._branch_block_dict = branch_block_dict
        # dictionary of {str: List of blocks}
        self._branch_block_dict_python = {i: x for i, x in branch_block_dict.items()}
    
    def _get_Sailfish_blocks(self) -> Dict[str, _Sailfish.Block]:
        return self._branch_block_dict
    
    def get_branches_str(self) -> str:
        return {i: self._branch_block_dict[i].print_tree() for i in list(self._branch_block_dict.keys())}
    
    def get_specific_branch(self, branch: str) -> str:
        if not branch in self._branch_block_dict_python:
            raise ValueError(f"branch not in the _branch_block, aviable branches are: {list(self._branch_block_dict_python.keys())}")
        return self._branch_block_dict[branch].print_tree()
    
    def print_branches(self) -> str:
        for i in list(self._branch_block_dict.keys()):
            print(f"branch = {i}")
            print(self._branch_block_dict[i].print_tree())
    
    def block_list(self)  -> List:
        if not branch in self._branch_block_dict_python:
            raise ValueError(f"branch not in the _branch_block, aviable branches are: {list(self._branch_block_dict_python.keys())}")
        return self._branch_block_dict_python[branch]

class Simulator:
    """Simulate MSAs based on SimProtocol"""
    
    def __init__(
        self, 
        simProtocol: Optional[SimProtocol] = None,
        simulation_type: Optional[SIMULATION_TYPE] = None
    ):
        if not simProtocol:
            warnings.warn("initalized a simulator without simProtocol -> using a default protocol with Tree = '(A:0.01,B:0.5,C:0.03);' and root length of 100")
            # default simulation values
            possion = PoissonDistribution(10, 100)
            simProtocol = SimProtocol(tree="(A:0.01,B:0.5,C:0.03);")
            simProtocol.set_insertion_length_distributions(possion)
            simProtocol.set_deletion_length_distributions(possion)
            simProtocol.set_insertion_rates(0.05)
            simProtocol.set_deletion_rates(0.05)
            simProtocol.set_sequence_size(100)
            simProtocol.set_min_sequence_size(1)

        # verify sim_protocol
        if self._verify_sim_protocol(simProtocol):
            self._simProtocol = simProtocol
            if simulation_type == SIMULATION_TYPE.PROTEIN:
                self._simulator = _Sailfish.AminoSimulator(self._simProtocol._sim)
            else:
                self._simulator = _Sailfish.NucleotideSimulator(self._simProtocol._sim)
        else:
            raise ValueError("failed to verify simProtocol")
        
        if not simulation_type:
            warnings.warn("simulation type not provided -> running indel only simulation")
            simulation_type = SIMULATION_TYPE.NOSUBS
        
        if simulation_type == SIMULATION_TYPE.PROTEIN:
            self._alphabet = _Sailfish.alphabetCode.AMINOACID
        elif simulation_type == SIMULATION_TYPE.DNA:
            self._alphabet = _Sailfish.alphabetCode.NUCLEOTIDE
        elif simulation_type == SIMULATION_TYPE.NOSUBS:
            self._alphabet = _Sailfish.alphabetCode.NULLCODE
        else:
            raise ValueError(f"unknown simulation type, please provde one of the following: {[e.name for e in SIMULATION_TYPE]}")
        
        self._simulation_type = simulation_type
        self._is_sub_model_init = False
    
    def _verify_sim_protocol(self, simProtocol) -> bool:
        if not simProtocol.get_tree():
            raise ValueError("protocol miss tree, please provide when initalizing the simProtocol")
        if not simProtocol.get_sequence_size() or simProtocol.get_sequence_size() == 0:
            raise ValueError("protocol miss root length, please provide -> simProtocol.set_sequence_size(int)")
        if not simProtocol.get_insertion_length_distribution(0):
            raise ValueError("protocol miss insertion length distribution, please provide -> simProtocol.set_insertion_length_distributions(float)")
        if not simProtocol.get_deletion_length_distribution(0):
            raise ValueError("protocol miss deletion length distribution, please provide -> simProtocol.set_deletion_length_distributions(float)")
        if simProtocol.get_insertion_rate(0) < 0:
            raise ValueError(f"please provide a non zero value for insertion rate, provided value of: {simProtocol.get_insertion_rate(0)} -> simProtocol.set_insertion_rate(float)")
        if simProtocol.get_deletion_rate(0) < 0:
            raise ValueError(f"please provide a non zero value for deletion rate, provided value of: {simProtocol.get_deletion_rate(0)} -> simProtocol.set_deletion_rate(float)")
        return True
    
    def reset_sim(self):
        # TODO, complete
        pass

    def _create_site_rate_model(
        self,
        gamma_alpha: float = 1.0,
        gamma_categories: int = 1,
        invariant_proportion: float = 0.0,
        site_rate_correlation: float = 0.0
    ) -> tuple[List[float], List[float], List[List[float]]]:
        """
        Create rate categories and probabilities for the site rate model.
        
        Args:
            gamma_alpha: Alpha parameter for gamma distribution
            gamma_categories: Number of gamma rate categories
            invariant_proportion: Proportion of invariant sites (0 to <1)
            site_rate_correlation: Correlation between adjacent sites (0 to <1)
                This is the ρ parameter for bivariate normal correlation.
                The realized discrete gamma correlation ρ_dG will be somewhat different.
        
        Returns:
            Tuple of (rates, probabilities, transition_matrix)
            - rates: List of rate values for each category
            - probabilities: Stationary probabilities for each category
            - transition_matrix: K×K matrix for correlated rates (empty list if no correlation)
        """
        if invariant_proportion < 0.0 or invariant_proportion >= 1.0:
            raise ValueError(f"invariant_proportion must be in [0, 1), received: {invariant_proportion}")
        
        if site_rate_correlation < 0.0 or site_rate_correlation >= 1.0:
            raise ValueError(f"site_rate_correlation must be in [0, 1), received: {site_rate_correlation}")
        
        # Validate correlation requires multiple categories
        if site_rate_correlation > 0.0 and gamma_categories == 1:
            warnings.warn(
                "site_rate_correlation > 0 requires gamma_categories > 1. "
                "Setting site_rate_correlation to 0.0"
            )
            site_rate_correlation = 0.0
        
        # Create gamma distribution
        gamma_dist = _Sailfish.GammaDistribution(gamma_alpha, gamma_categories)
        rates = list(gamma_dist.getAllRates())
        probs = list(gamma_dist.getAllRatesProb())
        
        # Add invariant sites category if requested
        if invariant_proportion > 0.0:
            # Scale existing probabilities
            scale_factor = 1.0 - invariant_proportion
            probs = [p * scale_factor for p in probs]
            
            # Add invariant category at the beginning
            rates.insert(0, 0.0)
            probs.insert(0, invariant_proportion)
        
        # Build transition matrix for correlated rates
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
                        rho=site_rate_correlation
                    )
                    
                except ImportError:
                    warnings.warn(
                        "site_rate_correlation > 0 requires scipy. "
                        "Install with: pip install scipy or pip install 'msasim[correlation]'. "
                        "Ignoring correlation parameter."
                    )
        return rates, probs, transition_matrix
    
    def _init_sub_model(self) -> None:
        self._model_factory = _Sailfish.modelFactory(self._simProtocol._get_Sailfish_tree())
        self._model_factory.set_alphabet(self._alphabet)
        if self._simulation_type == SIMULATION_TYPE.PROTEIN:
            warnings.warn("replacement matrix not provided -> running with default parameters: WAG model")
            self._model_factory.set_replacement_model(_Sailfish.modelCode.WAG)
        else:
            warnings.warn("replacement matrix not provided -> running with default parameters: JC model")
            self._model_factory.set_replacement_model(_Sailfish.modelCode.NUCJC)

    
        rates, probs = self._create_site_rate_model()
        self._model_factory.setSiteRateModel(rates, probs)    



        self._simulator.init_substitution_sim(self._model_factory)
        self._is_sub_model_init = True
    
    def set_replacement_model(
            self,
            model: _Sailfish.modelCode,
            amino_model_file: pathlib.Path = None,
            model_parameters: List = None,
            gamma_parameters_alpha : float = 1.0,
            gamma_parameters_categories: int = 1,
            invariant_sites_proportion: float = 0.0,
            site_rate_correlation: float = 0.0,
        ) -> None:
        if not model:
            raise ValueError(f"please provide a substitution model from the the following list: {_Sailfish.modelCode}")
        if int(gamma_parameters_categories) != gamma_parameters_categories:
            raise ValueError(f"gamma_parameters_catergories has to be a positive int value: received value of {gamma_parameters_categories}")
        self._model_factory = _Sailfish.modelFactory(self._simProtocol._get_Sailfish_tree())

        self._model_factory.set_alphabet(self._alphabet)
        if self._simulation_type == SIMULATION_TYPE.PROTEIN:
            if model_parameters:
                raise ValueError(f"no model parameters are used in protein, recevied value of: {model_parameters}")
            self._model_factory.set_replacement_model(model)
            if model == MODEL_CODES.CUSTOM and amino_model_file:
                self._model_factory.set_amino_replacement_model_file(str(amino_model_file))
        else:
            if model == MODEL_CODES.NUCJC and model_parameters:
                raise ValueError(f"no model parameters in JC model, recevied value of: {model_parameters}")
            self._model_factory.set_replacement_model(model)
            if model == MODEL_CODES.NUCJC and not model_parameters:
                pass
            elif not model_parameters:
                raise ValueError("please provide a model parameters")
            else:
                self._model_factory.set_model_parameters(model_parameters)


        rates, probs, transition_matrix = self._create_site_rate_model(
            gamma_alpha=gamma_parameters_alpha,
            gamma_categories=gamma_parameters_categories,
            invariant_proportion=invariant_sites_proportion,
            site_rate_correlation=site_rate_correlation
        )
        self._model_factory.setSiteRateModel(rates, probs, transition_matrix)

        self._simulator.init_substitution_sim(self._model_factory)

        self._is_sub_model_init = True
    
    def gen_indels(self) -> BlockTreePython:
        return BlockTreePython(self._simulator.gen_indels())
    
    def get_sequences_to_save(self) -> List[bool]:
        return self._simulator.get_saved_nodes_mask()
    
    def save_root_sequence(self):
        self._simulator.save_root_sequence()
    
    def save_all_nodes_sequences(self):
        self._simulator.save_all_nodes_sequences()

    def gen_substitutions(self, length: int):
        if not self._is_sub_model_init:
            self._init_sub_model()
        return self._simulator.gen_substitutions(length)
    
    # @profile
    def simulate(self, times: int = 1) -> List[Msa]:
        Msas = []
        for _ in range(times):
            if self._simProtocol._is_insertion_rate_zero and self._simProtocol._is_deletion_rate_zero:
                msa = Msa(sum(self.get_sequences_to_save()),
                          self._simProtocol.get_sequence_size(),
                          self.get_sequences_to_save())
            else:
                blocktree = self.gen_indels()
                msa = Msa(blocktree._get_Sailfish_blocks(),
                          self._simProtocol._get_root(),
                          self.get_sequences_to_save())

            # sim.init_substitution_sim(mFac)
            if self._simulation_type != SIMULATION_TYPE.NOSUBS:
                substitutions = self.gen_substitutions(msa.get_length())
                msa.fill_substitutions(substitutions)

            Msas.append(msa)
        return Msas
    
    def simulate_low_memory(self, output_file_path: pathlib.Path) -> Msa:
        if self._simProtocol._is_insertion_rate_zero and self._simProtocol._is_deletion_rate_zero:
            msa_length = self._simProtocol.get_sequence_size()
        else:
            blocktree = self.gen_indels()
            msa = Msa(blocktree._get_Sailfish_blocks(),
                        self._simProtocol._get_root(),
                        self.get_sequences_to_save())
            msa_length = msa.get_length()
            self._simulator.set_aligned_sequence_map(msa._msa)

        # sim.init_substitution_sim(mFac)
        if self._simulation_type != SIMULATION_TYPE.NOSUBS:
            self._simulator.gen_substitutions_to_file(msa_length, 
                                                      str(output_file_path))
        else:
            msa.write_msa(str(output_file_path))
    
    def __call__(self) -> Msa:
        return self.simulate(1)[0]
    
    def save_rates(self, is_save: bool) -> None:
        self._simulator.save_site_rates(is_save)
    
    def get_rates(self) -> List[float]:
        return self._simulator.get_site_rates()