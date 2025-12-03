"""Simulation protocol configuration"""

import _Sailfish
from typing import List, Optional, Dict
from .tree import Tree
from .distributions import Distribution, ZipfDistribution

class SimProtocol:
    """Configuration for MSA simulation"""
    
    def __init__(
        self, 
        tree = None,
        root_seq_size: int = 100,
        deletion_rate: float = 0.0,
        insertion_rate: float = 0.0,
        deletion_dist: Distribution = None,
        insertion_dist: Distribution = None,
        minimum_seq_size: int = 100,
        seed: int = 0,
    ):
        # Handle defaults
        if deletion_dist is None:
            deletion_dist = ZipfDistribution(1.7, 50)
        if insertion_dist is None:
            insertion_dist = ZipfDistribution(1.7, 50)
            
        # Parse tree
        if isinstance(tree, Tree):
            self._tree = tree
        elif isinstance(tree, str):
            self._tree = Tree(tree)
        else:
            raise ValueError("tree must be Tree object or newick string/path")
        
        self._num_branches = self._tree.get_num_nodes() - 1
        self._sim = _Sailfish.SimProtocol(self._tree._get_Sailfish_tree())
        self.set_seed(seed)
        self.set_sequence_size(root_seq_size)
        self._is_deletion_rate_zero = not deletion_rate
        self._is_insertion_rate_zero = not insertion_rate
        self.set_deletion_rates(deletion_rate=deletion_rate)
        self.set_insertion_rates(insertion_rate=insertion_rate)
        self.set_deletion_length_distributions(deletion_dist=deletion_dist)
        self.set_insertion_length_distributions(insertion_dist=insertion_dist)
        self.set_min_sequence_size(min_sequence_size=minimum_seq_size)

    def get_tree(self) -> Tree:
        return self._tree
    
    def _get_Sailfish_tree(self) -> _Sailfish.Tree:
        return self._tree._get_Sailfish_tree()
    
    def _get_root(self):
        return self._tree._get_Sailfish_tree().root
    
    def get_num_branches(self) -> int:
        return self._num_branches
    
    def set_seed(self, seed: int) -> None:
        self._seed = seed
        self._sim.set_seed(seed)
    
    def get_seed(self) -> int:
        return self._seed
    
    def set_sequence_size(self, sequence_size: int) -> None:
        self._sim.set_sequence_size(sequence_size)
        self._root_seq_size = sequence_size
    
    def get_sequence_size(self) -> int:
        return self._root_seq_size
    
    def set_min_sequence_size(self, min_sequence_size: int) -> None:
        self._sim.set_minimum_sequence_size(min_sequence_size)
        self._min_seq_size = min_sequence_size

    
    def set_insertion_rates(self, insertion_rate: Optional[float] = None, insertion_rates: Optional[List[float]] = None) -> None:
        if insertion_rate is not None:
            self.insertion_rates = [insertion_rate] * self._num_branches
            if insertion_rate:
                self._is_insertion_rate_zero = False
        elif insertion_rates:
            if not len(insertion_rates) == self._num_branches:
                raise ValueError(f"The length of the insertaion rates should be equal to the number of branches in the tree. The insertion_rates length is {len(insertion_rates)} and the number of branches is {self._num_branches}. You can pass a single value as insertion_rate which will be used for all branches.")
            self.insertion_rates = insertion_rates
            for insertion_rate in insertion_rates:
                if insertion_rate:
                    self._is_insertion_rate_zero = False
        else:
            raise ValueError("please provide one of the following: insertion_rate (a single value used for all branches), or a insertion_rates (a list of values, each corresponding to a different branch)")
        
        self._sim.set_insertion_rates(self.insertion_rates)
    
    def get_insertion_rate(self, branch_num: int) -> float:
        if branch_num >= self._num_branches:
            raise ValueError(f"The branch number should be between 0 to {self._num_branches} (not included). Received value of {branch_num}")
        return self._sim.get_insertion_rate(branch_num)
    
    def get_all_insertion_rates(self) -> Dict:
        return {i: self.get_insertion_rate(i) for i in range(self._num_branches)}
    
    def set_deletion_rates(self, deletion_rate: Optional[float] = None, deletion_rates: Optional[List[float]] = None) -> None:
        if deletion_rate is not None:
            self.deletion_rates = [deletion_rate] * self._num_branches
            if deletion_rate:
                self._is_deletion_rate_zero = False
        elif deletion_rates:
            if not len(deletion_rates) == self._num_branches:
                raise ValueError(f"The length of the deletion rates should be equal to the number of branches in the tree. The deletion_rates length is {len(deletion_rates)} and the number of branches is {self._num_branches}. You can pass a single value as deletion_rate which will be used for all branches.")
            self.deletion_rates = deletion_rates
            for deletion_rate in deletion_rates:
                if deletion_rate:
                    self._is_deletion_rate_zero = False
        else:
            raise ValueError("please provide one of the following: deletion_rate (a single value used for all branches), or a deletion_rates (a list of values, each corresponding to a different branch)")
        
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
            raise ValueError("please provide one of the following: deletion_rate (a single value used for all branches), or a deletion_rates (a list of values, each corresponding to a different branch)")
        
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
            raise ValueError("please provide one of the following: deletion_rate (a single value used for all branches), or a deletion_rates (a list of values, each corresponding to a different branch)")
        
        self._sim.set_deletion_length_distributions([dist._get_Sailfish_dist() for dist in self.deletion_dists])
    
    def get_deletion_length_distribution(self, branch_num: int) -> Distribution:
        if branch_num >= self._num_branches:
            raise ValueError(f"The branch number should be between 0 to {self._num_branches} (not included). Received value of {branch_num}")
        return self.deletion_dists[branch_num]
    
    def get_all_deletion_length_distribution(self) -> Dict:
        return {i: self.get_deletion_length_distribution(i) for i in range(self._num_branches)}
