import _Sailfish
import os
import warnings

from typing import List, Optional, Dict
from re import split
from enum import Enum

import numpy as np
from scipy.stats import geom, poisson, zipf

print(_Sailfish.Tree)
print(_Sailfish.Simulator)
print(_Sailfish.SimProtocol)
print(_Sailfish.Msa)
print(_Sailfish.DiscreteDistribution)
print(_Sailfish.modelFactory)
print(_Sailfish.alphabetCode)
print(_Sailfish.modelCode)
 

class SIMULATION_TYPE(Enum):
    DNA = 1
    PROTEIN = 2

def is_newick(tree: str):
    # from: https://github.com/ila/Newick-validator/blob/master/Newick_Validator.py
    # dividing the string into tokens, to check them singularly
    tokens = split(r'([A-Za-z]+[^A-Za-z,)]+[A-Za-z]+|[0-9.]*[A-Za-z]+[0-9.]+|[0-9.]+\s+[0-9.]+|[0-9.]+|[A-za-z]+|\(|\)|;|:|,)', tree)

    # removing spaces and empty strings (spaces within labels are still present)
    parsed_tokens = list(filter(lambda x: not (x.isspace() or not x), tokens))

    # checking whether the tree ends with ;
    if parsed_tokens[-1] != ';':
        raise ValueError(f"Tree without ; at the end. Tree received: {tree}")
        return False
    return True

# TODO, I think should be deleted
class Block:
    '''
    A single block of event. 
    Used to add insertions or deletions.
    '''
    def __init__(self, num1: int, num2: int):
        self.block = _Sailfish.Block(num1, num2)

class BlockTree:
    '''
    Used to contain the events on a multiple branches (entire tree).
    '''
    def __init__(self, root_length: int):
        self.blockTree = _Sailfish.BlockTree(root_length)
    
    def print_tree(self) -> str:
        return self.blockTree.print_tree()
    
    def block_list(self)  -> List:
        return self.blockTree.block_list()

# TODO delete one of this (I think the above if not used)
class BlockTreePython:
    '''
    Used to contain the events on a multiple branches (entire tree).
    '''
    def __init__(self, branch_block_dict: Dict[str, _Sailfish.Block]):
        self._branch_block_dict = branch_block_dict
        # dictionary of {str: List of blocks}
        self._branch_block_dict_python = {i: x.block_list() for i, x in branch_block_dict.items()}
    
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

class Distribution:
    def set_dist(self, dist):
        # sum should be "around" 1
        epsilon = 10e-6
        if abs(sum(dist)-1) > epsilon:
            raise ValueError(f"Sum of the distribution should be 1 for a valid probability distribution. Input received is: {dist}, sum is {sum(dist)}")
        for x in dist:
            if x < 0 or x > 1:
                raise ValueError(f"Each value of the probabilities should be between 0 to 1. Received a value of {x}")
        self._dist = _Sailfish.DiscreteDistribution(dist)
    
    def draw_sample(self) -> int:
        return self._dist.draw_sample()
    
    def set_seed(self, seed: int) -> None:
        return self._dist.set_seed(seed)
    
    def get_table(self) -> List:
        return self._dist.get_table()
    
    def _get_Sailfish_dist(self) -> _Sailfish.DiscreteDistribution:
        return self._dist

class CustomDistribution(Distribution):
    '''
    Provide a custom discrete distribution to the model.
    '''
    def __init__(self, dist: List[float]):
        self.set_dist(dist)

class Tree:
    '''
    The tree class for the simulator
    '''
    def __init__(self, input_str: str):
        is_from_file = False
        if os.path.isfile(input_str):
            is_from_file = True
            tree_str = open(input_str, 'r').read()
        else:
            tree_str = input_str
        if not is_newick(tree_str):
            if is_from_file:
                raise ValueError(f"Failed to read tree from file. File path: {input_str}, content: {tree_str}")
            else:
                raise ValueError(f"Failed construct tree from string. String received: {tree_str}")
        self._tree = _Sailfish.Tree(input_str, is_from_file)
        self._tree_str = tree_str
    
    def get_num_nodes(self) -> int:
        return self._tree.num_nodes

    def _get_Sailfish_tree(self) -> _Sailfish.Tree:
        return self._tree
    
    def __repr__(self) -> str:
        return f"{self._tree_str}"

class SimProtocol:
    '''
    The simulator protocol, sets the different distribution, tree and root length.
    '''
    def __init__(self, tree = None):
        if isinstance(tree, Tree):
            self._tree = tree
        elif isinstance(tree, str):
            self._tree = Tree(tree) 
        else:
            raise ValueError(f"please provide one of the following: (1) a newick format of a tree; (2) a path to a file containing a tree; (3) or a tree created by the Tree class")
        
        self._num_branches = self._tree.get_num_nodes() - 1
        self._sim = _Sailfish.SimProtocol(self._tree._get_Sailfish_tree())
        self._root_seq_size = 0
    
    def get_tree(self) -> Tree:
        return self._tree
    
    def _get_Sailfish_tree(self) -> _Sailfish.Tree:
        return self._tree._get_Sailfish_tree()
    
    def _get_root(self):
        return self._tree._get_Sailfish_tree().root
    
    def get_num_branches(self) -> int:
        return self._num_branches
    
    def set_seed(self, seed: int) -> None:
        self._sim.set_seed(seed)
    
    def get_seed(self) -> int:
        return sim.get_seed()
    
    def set_sequence_size(self, sequence_size: int) -> None:
        self._sim.set_sequence_size(sequence_size)
        self._root_seq_size = sequence_size
    
    def get_sequence_size(self) -> int:
        return self._root_seq_size
    
    def set_insertion_rates(self, insertion_rate: Optional[float] = None, insertion_rates: Optional[List[float]] = None) -> None:
        if insertion_rate:
            self.insertion_rates = [insertion_rate] * self._num_branches
        elif insertion_rates:
            if not len(insertion_rates) == self._num_branches:
                raise ValueError(f"The length of the insertaion rates should be equal to the number of branches in the tree. The insertion_rates length is {len(insertion_rates)} and the number of branches is {self._num_branches}. You can pass a single value as insertion_rate which will be used for all branches.")
            self.insertion_rates = insertion_rates
        else:
            raise ValueError(f"please provide one of the following: insertion_rate (a single value used for all branches), or a insertion_rates (a list of values, each corresponding to a different branch)")
        
        self._sim.set_insertion_rates(self.insertion_rates)
    
    def get_insertion_rate(self, branch_num: int) -> float:
        if branch_num >= self._num_branches:
            raise ValueError(f"The branch number should be between 0 to {self._num_branches} (not included). Received value of {branch_num}")
        return self._sim.get_insertion_rate(branch_num)
    
    def get_all_insertion_rates(self) -> Dict:
        return {i: self.get_insertion_rate(i) for i in range(self._num_branches)}
    
    def set_deletion_rates(self, deletion_rate: Optional[float] = None, deletion_rates: Optional[List[float]] = None) -> None:
        if deletion_rate:
            self.deletion_rates = [deletion_rate] * self._num_branches
        elif deletion_rates:
            if not len(deletion_rates) == self._num_branches:
                raise ValueError(f"The length of the deletion rates should be equal to the number of branches in the tree. The deletion_rates length is {len(deletion_rates)} and the number of branches is {self._num_branches}. You can pass a single value as deletion_rate which will be used for all branches.")
            self.deletion_rates = deletion_rates
        else:
            raise ValueError(f"please provide one of the following: deletion_rate (a single value used for all branches), or a deletion_rates (a list of values, each corresponding to a different branch)")
        
        self._sim.set_deletion_rates(self.deletion_rates)
    
    def get_deletion_rate(self, branch_num: int) -> float:
        if branch_num >= self._num_branches:
            raise ValueError(f"The branch number should be between 0 to {self._num_branches} (not included). Received value of {branch_num}")
        return self._sim.get_deletion_rate(branch_num)
    
    def get_all_deletion_rates(self) -> Dict:
        return {i: self.get_deletion_rate(i) for i in range(self._num_branches)}
    
    def set_insertion_length_distributions(self, insertion_dist: Optional[Distribution] = None, insertion_dists: Optional[List[Distribution]] = None) -> None:
        if insertion_dist:
            self.insertion_dists = [insertion_dist] * self._num_branches
        elif insertion_dists:
            if not len(insertion_dists) == self._num_branches:
                raise ValueError(f"The length of the insertion dists should be equal to the number of branches in the tree. The insertion_dists length is {len(insertion_dists)} and the number of branches is {self._num_branches}. You can pass a single value as insertion_dist which will be used for all branches.")
            self.insertion_dists = insertion_dists
        else:
            raise ValueError(f"please provide one of the following: deletion_rate (a single value used for all branches), or a deletion_rates (a list of values, each corresponding to a different branch)")
        
        self._sim.set_insertion_length_distributions([dist._get_Sailfish_dist() for dist in self.insertion_dists])
    
    def get_insertion_length_distribution(self, branch_num: int) -> Distribution:
        if branch_num >= self._num_branches:
            raise ValueError(f"The branch number should be between 0 to {self._num_branches} (not included). Received value of {branch_num}")
        return self.insertion_dists[branch_num]
    
    def get_all_insertion_length_distribution(self) -> Dict:
        return {i: self.get_insertion_length_distribution(i) for i in range(self._num_branches)}
    
    def set_deletion_length_distributions(self, deletion_dist: Optional[Distribution] = None, deletion_dists: Optional[List[Distribution]] = None) -> None:
        if deletion_dist:
            self.deletion_dists = [deletion_dist] * self._num_branches
        elif deletion_dists:
            if not len(deletion_dists) == self._num_branches:
                raise ValueError(f"The length of the deletion dists should be equal to the number of branches in the tree. The deletion_dists length is {len(deletion_dists)} and the number of branches is {self._num_branches}. You can pass a single value as deletion_dist which will be used for all branches.")
            self.deletion_dists = deletion_dists
        else:
            raise ValueError(f"please provide one of the following: deletion_rate (a single value used for all branches), or a deletion_rates (a list of values, each corresponding to a different branch)")
        
        self._sim.set_deletion_length_distributions([dist._get_Sailfish_dist() for dist in self.deletion_dists])
    
    def get_deletion_length_distribution(self, branch_num: int) -> Distribution:
        if branch_num >= self._num_branches:
            raise ValueError(f"The branch number should be between 0 to {self._num_branches} (not included). Received value of {branch_num}")
        return self.deletion_dists[branch_num]
    
    def get_all_deletion_length_distribution(self) -> Dict:
        return {i: self.get_deletion_length_distribution(i) for i in range(self._num_branches)}

class Msa:
    '''
    The MSA class from the simulator
    '''
    def __init__(self, species_dict: Dict[str, BlockTree], root_node):
        self._msa = _Sailfish.Msa(species_dict, root_node)
    
    def generate_msas(self, node):
        self._msa.generate_msas(node)
    
    def get_length(self) -> int:
        return self._msa.length()
    
    def get_num_sequences(self) -> int:
        return self._msa.num_sequences()
    
    def fill_substitutions(self, sequenceContainer) -> None:
        self._msa.fill_substitutions(sequenceContainer)
    
    def print_msa(self) -> str:
        return self._msa.print_msa()
    
    def print_indels(self) -> str:
        return self._msa.print_indels()
    
    def get_msa(self) -> str:
        return self._msa.get_msa()
    
    #def __repr__(self) -> str:
    #    return f"{self.get_msa()}"

class Simulator:
    '''
    Simulate MSAs based on SimProtocol
    '''
    def __init__(self, simProtocol: Optional[SimProtocol] = None, simulation_type: Optional[SIMULATION_TYPE] = None):
        if not simProtocol:
            warnings.warn(f"initalized a simulator without simProtocol -> using a default protocol with Tree = '(A:0.01,B:0.5,C:0.03);' and root length of 100")
            # default simulation values
            possion = PoissonDistribution(10, 100)
            simProtocol = SimProtocol(tree="(A:0.01,B:0.5,C:0.03);")
            simProtocol.set_insertion_length_distributions(possion)
            simProtocol.set_deletion_length_distributions(possion)
            simProtocol.set_insertion_rates(0.05)
            simProtocol.set_deletion_rates(0.05)
            simProtocol.set_sequence_size(100)
        # verify sim_protocol
        if self._verify_sim_protocol(simProtocol):
            self._simProtocol = simProtocol
            self._simulator = _Sailfish.Simulator(self._simProtocol._sim)
        else:
            raise ValueError(f"failed to verify simProtocol")
        
        if not simulation_type:
            warnings.warn(f"simulation type not provided -> running protein simulation")
            simulation_type = SIMULATION_TYPE.PROTEIN
        
        if simulation_type == SIMULATION_TYPE.PROTEIN:
            self._alphabet = _Sailfish.alphabetCode.AMINOACID
        elif simulation_type == SIMULATION_TYPE.DNA:
            self._alphabet = _Sailfish.alphabetCode.NUCLEOTIDE
        else:
            raise ValueError(f"unknown simulation type, please provde one of the following: {[e.name for e in SIMULATION_TYPE]}")
        
        self._simulation_type = simulation_type
        self._is_sub_model_init = False
    
    def _verify_sim_protocol(self, simProtocol) -> bool:
        if not simProtocol.get_tree():
            raise ValueError(f"protocol miss tree, please provide when initalizing the simProtocol")
        if not simProtocol.get_sequence_size() or simProtocol.get_sequence_size() == 0:
            raise ValueError(f"protocol miss root length, please provide -> simProtocol.set_sequence_size(int)")
        if not simProtocol.get_insertion_length_distribution(0):
            raise ValueError(f"protocol miss insertion length distribution, please provide -> simProtocol.set_insertion_length_distributions(float)")
        if not simProtocol.get_deletion_length_distribution(0):
            raise ValueError(f"protocol miss deletion length distribution, please provide -> simProtocol.set_deletion_length_distributions(float)")
        if not simProtocol.get_insertion_rate(0):
            raise ValueError(f"protocol miss insertion rate, please provide -> simProtocol.set_insertion_rate(float)")
        if not simProtocol.get_deletion_rate(0):
            raise ValueError(f"protocol miss deletion rate, please provide -> simProtocol.set_deletion_rate(float)")
        return True
    
    def reset_sim(self):
        # TODO, complete
        pass
    
    def _init_sub_model(self) -> None:
        self._model_factory = _Sailfish.modelFactory(self._simProtocol._get_Sailfish_tree())
        self._model_factory.set_alphabet(self._alphabet)
        if self._simulation_type == SIMULATION_TYPE.PROTEIN:
            warnings.warn(f"replacement matrix not provided -> running with default parameters: WAG model code with gamma parameters of (0.5, 4)")
            self._model_factory.set_replacement_model(_Sailfish.modelCode.WAG)
            self._model_factory.set_gamma_parameters(0.5, 4)
        else:
            # TODO complete
            pass
        self._is_sub_model_init = True
    
    def set_replacement_model(model: _Sailfish.modelCode, *args) -> None:
        if self._simulation_type == SIMULATION_TYPE.PROTEIN:
            self._model_factory.set_replacement_model(model)
            self._model_factory.set_gamma_parameters(*args)
        else:
            self._model_factory.set_replacement_model(model)
            self._model_factory.set_model_parameters(*args)
        self._is_sub_model_init = True
    
    def gen_indels(self) -> BlockTreePython:
        return BlockTreePython(self._simulator.gen_indels())
    
    def gen_substitutions(self, length: int):
        if not self._is_sub_model_init:
            self._init_sub_model()
        return self._simulator.gen_substitutions(self._model_factory, length)
    
    def simulate(self, times: int = 1) -> List[Msa]:
        Msas = []
        for idx in range(times):
            blocktree = self.gen_indels()
            msa = Msa(blocktree._get_Sailfish_blocks(), self._simProtocol._get_root())
            substitutions = self.gen_substitutions(msa.get_length())
            msa.fill_substitutions(substitutions)
            if times == 1:
                return msa
            Msas.append(msa)
        return Msas

class GeometricDistribution(Distribution):
    def __init__(self, p: float, truncation: int = 150):
        """
        Calculation of geoemtric moment
        inputs:
        p - p parameter of the geoemtric distribution
        truncation - (optional, by default 150) maximal value of the distribution
        """
        self.p = p
        self.truncation = truncation
        CDF = lambda x: geom.cdf(x, p)
        norm_factor = CDF(truncation) - CDF(0)

        probabilities = geom.pmf(np.arange(1, truncation+1), p)
        probabilities = probabilities / norm_factor

        self.set_dist(probabilities)

    def __repr__(self) -> str:
        return f"Geometric distribution: (p={self.p}, truncation{self.truncation})"

class PoissonDistribution(Distribution):
    def __init__(self, p: float, truncation: int = 150):
        """
        Calculation of geoemtric moment
        inputs:
        p - p parameter of the geoemtric distribution
        truncation - (optional, by default 150) maximal value of the distribution
        """
        self.p = p
        self.truncation = truncation
        CDF = lambda x: poisson.cdf(x, p)
        norm_factor = CDF(truncation) - CDF(0)

        probabilities = poisson.pmf(np.arange(1, truncation+1), p)
        probabilities = probabilities / norm_factor

        self.set_dist(probabilities)

    def __repr__(self) -> str:
        return f"Poisson distribution: (p={self.p}, truncation{self.truncation})"

class ZipfDistribution(Distribution):
    def __init__(self, p: float, truncation: int = 150):
        """
        Calculation of geoemtric moment
        inputs:
        p - p parameter of the geoemtric distribution
        truncation - (optional, by default 150) maximal value of the distribution
        """
        self.p = p
        self.truncation = truncation
        CDF = lambda x: zipf.cdf(x, p)
        norm_factor = CDF(truncation) - CDF(0)

        probabilities = zipf.pmf(np.arange(1, truncation+1), p)
        probabilities = probabilities / norm_factor

        self.set_dist(probabilities)
    
    def __repr__(self) -> str:
        return f"Zipf distribution: (p={self.p}, truncation{self.truncation})"