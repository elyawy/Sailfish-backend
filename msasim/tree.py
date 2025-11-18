"""Phylogenetic tree handling"""

import _Sailfish
import os
from re import split

def is_newick(tree: str) -> bool:
    """Validate newick format"""
    # from: https://github.com/ila/Newick-validator/blob/master/Newick_Validator.py
    # dividing the string into tokens, to check them singularly
    tokens = split(r'([A-Za-z]+[^A-Za-z,)]+[A-Za-z]+|[0-9.]*[A-Za-z]+[0-9.]+|[0-9.]+\s+[0-9.]+|[0-9.]+|[A-za-z]+|\(|\)|;|:|,)', tree)
    # removing spaces and empty strings (spaces within labels are still present)
    parsed_tokens = list(filter(lambda x: not (x.isspace() or not x), tokens))

    # checking whether the tree ends with ;
    if parsed_tokens[-1] != ';':
        raise ValueError(f"Tree must end with ';'. Got: {tree}")
    return True

class Tree:
    """Phylogenetic tree for simulation"""
    
    def __init__(self, input_str: str):
        is_from_file = os.path.isfile(input_str)
        
        if is_from_file:
            tree_str = open(input_str, 'r').read()
        else:
            tree_str = input_str
            
        if not is_newick(tree_str):
            source = "file" if is_from_file else "string"
            raise ValueError(f"Invalid newick from {source}: {tree_str}")
            
        self._tree = _Sailfish.Tree(input_str, is_from_file)
        self._tree_str = tree_str
    
    def get_num_nodes(self) -> int:
        return self._tree.num_nodes
    
    def get_num_leaves(self) -> int:
        return self._tree.root.num_leaves
    
    def _get_Sailfish_tree(self) -> _Sailfish.Tree:
        """Internal: Get C++ tree object"""
        return self._tree
    
    def __repr__(self) -> str:
        return self._tree_str