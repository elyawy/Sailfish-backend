import _Sailfish
import os

from typing import List, Optional, Dict
from re import split

print(_Sailfish.Tree)
print(_Sailfish.Simulator)
print(_Sailfish.SimProtocol)
print(_Sailfish.Msa)
print(_Sailfish.DiscreteDistribution)
print(_Sailfish.modelFactory)
print(_Sailfish.alphabetCode)
print(_Sailfish.modelCode)
 

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

class Block:
    '''
    A single block of event. 
    Used to add insertions or deletions.
    '''
    def __init__(self, num1: int, num2: int):
        self.block = _Sailfish.Block(num1, num2)

class BlockTree:
    '''
    Multiple blocks in a tree.
    Used to contain the events on a single branch.
    '''
    def __init__(self, root_length: int):
        self.blockTree = _Sailfish.BlockTree(root_length)
    
    def print_tree(self) -> str:
        return self.blockTree.print_tree()
    
    def block_list(self)  -> List:
        return self.blockTree.block_list()

class Distribution:
    def set_dist(self, dist):
        if sum(dist) != 1:
            raise ValueError(f"Sum of the distribution should be 1 for a valid probability distribution. Input received is: {dist}, sum is {sum(dist)}")
        for x in dist:
            if x < 0:
                raise ValueError(f"Each value of the probabilities should be between 0 to 1. Received a value of {x}")
        self.dist = _Sailfish.DiscreteDistribution(dist)
    
    def draw_sample(self) -> int:
        return self.dist.draw_sample()
    
    def set_seed(self, seed: int) -> None:
        return self.dist.set_seed(seed)
    
    def get_table(self) -> List:
        return self.dist.get_table()

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
        self.tree = _Sailfish.Tree(tree_str)

class SimProtocol:
    '''
    The simulator.
    '''
    def __init__(self, input_str: Optional[str] = None, tree: Optional[Tree] = None):
        if tree:
            self.tree = tree
        elif input_str:
            self.tree = Tree(input_str) 
        else:
            raise ValueError(f"please provide one of the following: input_str (newick format / pass to file containing a tree), or a tree (created by the Tree class)")
        
        self.num_branches = self.tree.num_nodes - 1
        self.sim = _Sailfish.SimulationProtocol(self.tree)
    
    def set_seed(self, seed: int) -> None:
        self.sim.set_seed(seed)
    
    def get_seed(self) -> int:
        return sim.get_seed()
    
    def set_sequence_size(self, sequence_size: int) -> None:
        self.sim.set_sequence_size(sequence_size)
    
    def get_sequence_size(self) -> int:
        return self.sim.get_sequence_size()
    
    def set_insertion_rates(self, insertion_rate: Optional[float] = None, insertion_rates: Optional[List[float]] = None) -> None:
        if insertion_rate:
            self.insertion_rates = [insertion_rate] * self.num_branches
        elif insertion_rates:
            if not len(insertion_rates) == self.num_branches:
                raise ValueError(f"The length of the insertaion rates should be equal to the number of branches in the tree. The insertion_rates length is {len(insertion_rates)} and the number of branches is {self.num_branches}. You can pass a single value as insertion_rate which will be used for all branches.")
            self.insertion_rates = insertion_rates
        else:
            raise ValueError(f"please provide one of the following: insertion_rate (a single value used for all branches), or a insertion_rates (a list of values, each corresponding to a different branch)")
        
        self.sim.set_insertion_rates(self.insertion_rates)
    
    def get_insertion_rate(self, branch_num: int) -> float:
        if branch_num >= self.num_branches:
            raise ValueError(f"The branch number should be between 0 to {self.num_branches} (not included). Received value of {branch_num}")
        return self.sim.get_insertion_rate(branch_num)
    
    def get_all_insertion_rates(self) -> Dict:
        return {i: self.get_insertion_rate(i) for i in range(self.num_branches)}
    
    def set_deletion_rates(self, deletion_rate: Optional[float] = None, deletion_rates: Optional[List[float]] = None) -> None:
        if deletion_rate:
            self.deletion_rates = [deletion_rate] * self.num_branches
        elif deletion_rates:
            if not len(deletion_rates) == self.num_branches:
                raise ValueError(f"The length of the deletion rates should be equal to the number of branches in the tree. The deletion_rates length is {len(deletion_rates)} and the number of branches is {self.num_branches}. You can pass a single value as deletion_rate which will be used for all branches.")
            self.deletion_rates = deletion_rates
        else:
            raise ValueError(f"please provide one of the following: deletion_rate (a single value used for all branches), or a deletion_rates (a list of values, each corresponding to a different branch)")
        
        self.sim.set_deletion_rates(self.deletion_rates)
    
    def get_deletion_rate(self, branch_num: int) -> float:
        if branch_num >= self.num_branches:
            raise ValueError(f"The branch number should be between 0 to {self.num_branches} (not included). Received value of {branch_num}")
        return self.sim.get_deletion_rate(branch_num)
    
    def get_all_deletion_rates(self) -> Dict:
        return {i: self.get_deletion_rate(i) for i in range(self.num_branches)}
    
    def set_insertion_length_distributions(self, insertion_dist: Optional[Distribution] = None, insertion_dists: Optional[List[Distribution]] = None) -> None:
        if insertion_dist:
            self.insertion_dists = [insertion_dist] * self.num_branches
        elif insertion_dists:
            if not len(insertion_dists) == self.num_branches:
                raise ValueError(f"The length of the insertion dists should be equal to the number of branches in the tree. The insertion_dists length is {len(insertion_dists)} and the number of branches is {self.num_branches}. You can pass a single value as insertion_dist which will be used for all branches.")
            self.insertion_dists = insertion_dists
        else:
            raise ValueError(f"please provide one of the following: deletion_rate (a single value used for all branches), or a deletion_rates (a list of values, each corresponding to a different branch)")
        
        self.sim.set_insertion_length_distributions(self.insertion_dists)
    
    def get_insertion_length_distribution(self, branch_num: int) -> Distribution:
        if branch_num >= self.num_branches:
            raise ValueError(f"The branch number should be between 0 to {self.num_branches} (not included). Received value of {branch_num}")
        return self.sim.get_insertion_length_distribution(branch_num)
    
    def get_all_insertion_length_distribution(self) -> Dict:
        return {i: self.get_insertion_length_distribution(i) for i in range(self.num_branches)}
    
    def set_deletion_length_distributions(self, deletion_dist: Optional[Distribution] = None, deletion_dists: Optional[List[Distribution]] = None) -> None:
        if deletion_dist:
            self.deletion_dists = [deletion_dist] * self.num_branches
        elif deletion_dists:
            if not len(deletion_dists) == self.num_branches:
                raise ValueError(f"The length of the deletion dists should be equal to the number of branches in the tree. The deletion_dists length is {len(deletion_dists)} and the number of branches is {self.num_branches}. You can pass a single value as deletion_dist which will be used for all branches.")
            self.deletion_dists = deletion_dists
        else:
            raise ValueError(f"please provide one of the following: deletion_rate (a single value used for all branches), or a deletion_rates (a list of values, each corresponding to a different branch)")
        
        self.sim.set_deletion_length_distributions(self.deletion_dists)
    
    def get_deletion_length_distribution(self, branch_num: int) -> Distribution:
        if branch_num >= self.num_branches:
            raise ValueError(f"The branch number should be between 0 to {self.num_branches} (not included). Received value of {branch_num}")
        return self.sim.get_deletion_length_distribution(branch_num)
    
    def get_all_deletion_length_distribution(self) -> Dict:
        return {i: self.get_deletion_length_distribution(i) for i in range(self.num_branches)}

class Msa:
    '''
    The MSA class from the simulator
    '''
    def __init__(self, species_dict: Dict[str, BlockTree], tree: Tree):
        self.msa = _Sailfish.Msa(species_dict, tree)
    
    def generate_msas(self, node):
        self.msa.generate_msas(node)
    
    def get_length(self) -> int:
        return self.msa.length()
    
    def get_num_sequences(self) -> int:
        return self.msa.num_sequences()
    
    def fill_substitutions(self, sequenceContainer) -> None:
        self.fill_substitutions(sequenceContainer)
    
    def print_msa(self) -> str:
        return self.msa.print_msa()
    
    def print_indels(self) -> str:
        return self.msa.print_indels()
    
    def get_msa(self) -> str:
        return self..msa.get_msa()