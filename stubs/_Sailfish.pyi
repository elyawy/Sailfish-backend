"""

        Sailfish simulator
        -----------------------

        .. currentmodule:: _Sailfish

        .. autosummary::
           :toctree: _generate

            DiscreteDistribution
            SimProtocol
            alphabetCode
            modelCode
            modelFactory
            Simulator
            Msa
            Tree
    
"""
from __future__ import annotations
import pybind11_stubgen.typing_ext
import typing
__all__ = ['AAJC', 'AMINOACID', 'Block', 'BlockTree', 'CPREV45', 'CUSTOM', 'DAYHOFF', 'Deletion', 'DiscreteDistribution', 'EHO_EXTENDED', 'EHO_HELIX', 'EHO_OTHER', 'EMPIRICODON', 'EX_BURIED', 'EX_EHO_BUR_EXT', 'EX_EHO_BUR_HEL', 'EX_EHO_BUR_OTH', 'EX_EHO_EXP_EXT', 'EX_EHO_EXP_HEL', 'EX_EHO_EXP_OTH', 'EX_EXPOSED', 'GTR', 'HIVB', 'HIVW', 'HKY', 'Insertion', 'JONES', 'LG', 'MTREV24', 'Msa', 'NUCJC', 'NUCLEOTIDE', 'NULLCODE', 'SimProtocol', 'Simulator', 'TAMURA92', 'Tree', 'WAG', 'alphabetCode', 'event', 'modelCode', 'modelFactory', 'node', 'sequenceContainer']
class Block:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: int, arg1: int) -> None:
        ...
class BlockTree:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: int) -> None:
        ...
    def block_list(self) -> list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]]:
        ...
    def print_tree(self) -> str:
        ...
class DiscreteDistribution:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    @staticmethod
    def set_seed(arg0: int) -> None:
        """
        Set seed for the random number generator
        """
    def __init__(self, arg0: list[float]) -> None:
        ...
    def draw_sample(self) -> int:
        """
        Draw a random sample according to the given distribution
        """
    def get_table(self) -> list[tuple[float, int]]:
        """
        Get Vose's alias table (useful for debugging)
        """
class Msa:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self, arg0: int, arg1: int, arg2: list[bool]) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: dict[int, tuple[list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]], int]], arg1: node, arg2: list[bool]) -> None:
        ...
    def fill_substitutions(self, arg0: sequenceContainer) -> None:
        ...
    def generate_msas(self: list[dict[int, tuple[list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]], int]]], arg0: node, arg1: list[bool]) -> list[Msa]:
        ...
    def get_msa(self) -> dict[int, list[int]]:
        ...
    def get_msa_string(self) -> str:
        ...
    def length(self) -> int:
        ...
    def num_sequences(self) -> int:
        ...
    def print_indels(self) -> None:
        ...
    def print_msa(self) -> None:
        ...
    def set_substitutions_folder(self, arg0: str) -> None:
        ...
    def write_msa(self, arg0: str) -> None:
        ...
    def write_msa_from_dir(self, arg0: str) -> None:
        ...
class SimProtocol:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: Tree) -> None:
        ...
    def get_deletion_length_distribution(self, arg0: int) -> DiscreteDistribution:
        ...
    def get_deletion_rate(self, arg0: int) -> float:
        ...
    def get_insertion_length_distribution(self, arg0: int) -> DiscreteDistribution:
        ...
    def get_insertion_rate(self, arg0: int) -> float:
        ...
    def get_minimum_sequence_size(self) -> int:
        ...
    def get_seed(self) -> int:
        ...
    def get_sequence_size(self) -> int:
        ...
    def set_deletion_length_distributions(self, arg0: list[DiscreteDistribution]) -> None:
        ...
    def set_deletion_rates(self, arg0: list[float]) -> None:
        ...
    def set_insertion_length_distributions(self, arg0: list[DiscreteDistribution]) -> None:
        ...
    def set_insertion_rates(self, arg0: list[float]) -> None:
        ...
    def set_minimum_sequence_size(self, arg0: int) -> None:
        ...
    def set_seed(self, arg0: int) -> None:
        ...
    def set_sequence_size(self, arg0: int) -> None:
        ...
class Simulator:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: SimProtocol) -> None:
        ...
    def gen_indels(self) -> dict[int, tuple[list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]], int]]:
        ...
    def gen_substitutions(self, arg0: int) -> sequenceContainer:
        ...
    def gen_substitutions_to_dir(self, arg0: int, arg1: str) -> None:
        ...
    def get_saved_nodes_mask(self) -> list[bool]:
        ...
    def get_site_rates(self) -> list[float]:
        ...
    def init_substitution_sim(self, arg0: modelFactory) -> None:
        ...
    def reset_sim(self, arg0: SimProtocol) -> None:
        ...
    def run_sim(self, arg0: int) -> list[dict[int, tuple[list[typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]], int]]]:
        ...
    def save_all_nodes_sequences(self) -> None:
        ...
    def save_root_sequence(self) -> None:
        ...
    def save_site_rates(self, arg0: bool) -> None:
        ...
class Tree:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: str, arg1: bool) -> None:
        """
        Create Phylogenetic tree object from newick formatted file
        """
    @property
    def num_nodes(self) -> int:
        ...
    @property
    def root(self) -> ...:
        ...
class alphabetCode:
    """
    Members:
    
      NULLCODE
    
      NUCLEOTIDE
    
      AMINOACID
    """
    AMINOACID: typing.ClassVar[alphabetCode]  # value = <alphabetCode.AMINOACID: 2>
    NUCLEOTIDE: typing.ClassVar[alphabetCode]  # value = <alphabetCode.NUCLEOTIDE: 1>
    NULLCODE: typing.ClassVar[alphabetCode]  # value = <alphabetCode.NULLCODE: 0>
    __members__: typing.ClassVar[dict[str, alphabetCode]]  # value = {'NULLCODE': <alphabetCode.NULLCODE: 0>, 'NUCLEOTIDE': <alphabetCode.NUCLEOTIDE: 1>, 'AMINOACID': <alphabetCode.AMINOACID: 2>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class event:
    """
    Members:
    
      Insertion
    
      Deletion
    """
    Deletion: typing.ClassVar[event]  # value = <event.Deletion: 1>
    Insertion: typing.ClassVar[event]  # value = <event.Insertion: 0>
    __members__: typing.ClassVar[dict[str, event]]  # value = {'Insertion': <event.Insertion: 0>, 'Deletion': <event.Deletion: 1>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class modelCode:
    """
    Members:
    
      NUCJC
    
      AAJC
    
      GTR
    
      HKY
    
      TAMURA92
    
      CPREV45
    
      DAYHOFF
    
      JONES
    
      MTREV24
    
      WAG
    
      HIVB
    
      HIVW
    
      LG
    
      EMPIRICODON
    
      EX_BURIED
    
      EX_EXPOSED
    
      EHO_EXTENDED
    
      EHO_HELIX
    
      EHO_OTHER
    
      EX_EHO_BUR_EXT
    
      EX_EHO_BUR_HEL
    
      EX_EHO_BUR_OTH
    
      EX_EHO_EXP_EXT
    
      EX_EHO_EXP_HEL
    
      EX_EHO_EXP_OTH
    
      CUSTOM
    """
    AAJC: typing.ClassVar[modelCode]  # value = <modelCode.AAJC: 1>
    CPREV45: typing.ClassVar[modelCode]  # value = <modelCode.CPREV45: 6>
    CUSTOM: typing.ClassVar[modelCode]  # value = <modelCode.CUSTOM: 26>
    DAYHOFF: typing.ClassVar[modelCode]  # value = <modelCode.DAYHOFF: 7>
    EHO_EXTENDED: typing.ClassVar[modelCode]  # value = <modelCode.EHO_EXTENDED: 17>
    EHO_HELIX: typing.ClassVar[modelCode]  # value = <modelCode.EHO_HELIX: 18>
    EHO_OTHER: typing.ClassVar[modelCode]  # value = <modelCode.EHO_OTHER: 19>
    EMPIRICODON: typing.ClassVar[modelCode]  # value = <modelCode.EMPIRICODON: 14>
    EX_BURIED: typing.ClassVar[modelCode]  # value = <modelCode.EX_BURIED: 15>
    EX_EHO_BUR_EXT: typing.ClassVar[modelCode]  # value = <modelCode.EX_EHO_BUR_EXT: 20>
    EX_EHO_BUR_HEL: typing.ClassVar[modelCode]  # value = <modelCode.EX_EHO_BUR_HEL: 21>
    EX_EHO_BUR_OTH: typing.ClassVar[modelCode]  # value = <modelCode.EX_EHO_BUR_OTH: 22>
    EX_EHO_EXP_EXT: typing.ClassVar[modelCode]  # value = <modelCode.EX_EHO_EXP_EXT: 23>
    EX_EHO_EXP_HEL: typing.ClassVar[modelCode]  # value = <modelCode.EX_EHO_EXP_HEL: 24>
    EX_EHO_EXP_OTH: typing.ClassVar[modelCode]  # value = <modelCode.EX_EHO_EXP_OTH: 25>
    EX_EXPOSED: typing.ClassVar[modelCode]  # value = <modelCode.EX_EXPOSED: 16>
    GTR: typing.ClassVar[modelCode]  # value = <modelCode.GTR: 2>
    HIVB: typing.ClassVar[modelCode]  # value = <modelCode.HIVB: 11>
    HIVW: typing.ClassVar[modelCode]  # value = <modelCode.HIVW: 12>
    HKY: typing.ClassVar[modelCode]  # value = <modelCode.HKY: 3>
    JONES: typing.ClassVar[modelCode]  # value = <modelCode.JONES: 8>
    LG: typing.ClassVar[modelCode]  # value = <modelCode.LG: 13>
    MTREV24: typing.ClassVar[modelCode]  # value = <modelCode.MTREV24: 9>
    NUCJC: typing.ClassVar[modelCode]  # value = <modelCode.NUCJC: 0>
    TAMURA92: typing.ClassVar[modelCode]  # value = <modelCode.TAMURA92: 4>
    WAG: typing.ClassVar[modelCode]  # value = <modelCode.WAG: 10>
    __members__: typing.ClassVar[dict[str, modelCode]]  # value = {'NUCJC': <modelCode.NUCJC: 0>, 'AAJC': <modelCode.AAJC: 1>, 'GTR': <modelCode.GTR: 2>, 'HKY': <modelCode.HKY: 3>, 'TAMURA92': <modelCode.TAMURA92: 4>, 'CPREV45': <modelCode.CPREV45: 6>, 'DAYHOFF': <modelCode.DAYHOFF: 7>, 'JONES': <modelCode.JONES: 8>, 'MTREV24': <modelCode.MTREV24: 9>, 'WAG': <modelCode.WAG: 10>, 'HIVB': <modelCode.HIVB: 11>, 'HIVW': <modelCode.HIVW: 12>, 'LG': <modelCode.LG: 13>, 'EMPIRICODON': <modelCode.EMPIRICODON: 14>, 'EX_BURIED': <modelCode.EX_BURIED: 15>, 'EX_EXPOSED': <modelCode.EX_EXPOSED: 16>, 'EHO_EXTENDED': <modelCode.EHO_EXTENDED: 17>, 'EHO_HELIX': <modelCode.EHO_HELIX: 18>, 'EHO_OTHER': <modelCode.EHO_OTHER: 19>, 'EX_EHO_BUR_EXT': <modelCode.EX_EHO_BUR_EXT: 20>, 'EX_EHO_BUR_HEL': <modelCode.EX_EHO_BUR_HEL: 21>, 'EX_EHO_BUR_OTH': <modelCode.EX_EHO_BUR_OTH: 22>, 'EX_EHO_EXP_EXT': <modelCode.EX_EHO_EXP_EXT: 23>, 'EX_EHO_EXP_HEL': <modelCode.EX_EHO_EXP_HEL: 24>, 'EX_EHO_EXP_OTH': <modelCode.EX_EHO_EXP_OTH: 25>, 'CUSTOM': <modelCode.CUSTOM: 26>}
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class modelFactory:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: Tree) -> None:
        ...
    def reset(self) -> None:
        ...
    def set_alphabet(self, arg0: alphabetCode) -> None:
        ...
    def set_amino_replacement_model_file(self, arg0: str) -> None:
        ...
    def set_gamma_parameters(self, arg0: float, arg1: int) -> None:
        ...
    def set_invariant_sites_proportion(self, arg0: float) -> None:
        ...
    def set_model_parameters(self, arg0: list[float]) -> None:
        ...
    def set_replacement_model(self, arg0: modelCode) -> None:
        ...
class node:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def distance_to_father(self) -> float:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def num_leaves(self) -> int:
        ...
    @property
    def sons(self) -> list[node]:
        ...
class sequenceContainer:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self) -> None:
        ...
AAJC: modelCode  # value = <modelCode.AAJC: 1>
AMINOACID: alphabetCode  # value = <alphabetCode.AMINOACID: 2>
CPREV45: modelCode  # value = <modelCode.CPREV45: 6>
CUSTOM: modelCode  # value = <modelCode.CUSTOM: 26>
DAYHOFF: modelCode  # value = <modelCode.DAYHOFF: 7>
Deletion: event  # value = <event.Deletion: 1>
EHO_EXTENDED: modelCode  # value = <modelCode.EHO_EXTENDED: 17>
EHO_HELIX: modelCode  # value = <modelCode.EHO_HELIX: 18>
EHO_OTHER: modelCode  # value = <modelCode.EHO_OTHER: 19>
EMPIRICODON: modelCode  # value = <modelCode.EMPIRICODON: 14>
EX_BURIED: modelCode  # value = <modelCode.EX_BURIED: 15>
EX_EHO_BUR_EXT: modelCode  # value = <modelCode.EX_EHO_BUR_EXT: 20>
EX_EHO_BUR_HEL: modelCode  # value = <modelCode.EX_EHO_BUR_HEL: 21>
EX_EHO_BUR_OTH: modelCode  # value = <modelCode.EX_EHO_BUR_OTH: 22>
EX_EHO_EXP_EXT: modelCode  # value = <modelCode.EX_EHO_EXP_EXT: 23>
EX_EHO_EXP_HEL: modelCode  # value = <modelCode.EX_EHO_EXP_HEL: 24>
EX_EHO_EXP_OTH: modelCode  # value = <modelCode.EX_EHO_EXP_OTH: 25>
EX_EXPOSED: modelCode  # value = <modelCode.EX_EXPOSED: 16>
GTR: modelCode  # value = <modelCode.GTR: 2>
HIVB: modelCode  # value = <modelCode.HIVB: 11>
HIVW: modelCode  # value = <modelCode.HIVW: 12>
HKY: modelCode  # value = <modelCode.HKY: 3>
Insertion: event  # value = <event.Insertion: 0>
JONES: modelCode  # value = <modelCode.JONES: 8>
LG: modelCode  # value = <modelCode.LG: 13>
MTREV24: modelCode  # value = <modelCode.MTREV24: 9>
NUCJC: modelCode  # value = <modelCode.NUCJC: 0>
NUCLEOTIDE: alphabetCode  # value = <alphabetCode.NUCLEOTIDE: 1>
NULLCODE: alphabetCode  # value = <alphabetCode.NULLCODE: 0>
TAMURA92: modelCode  # value = <modelCode.TAMURA92: 4>
WAG: modelCode  # value = <modelCode.WAG: 10>
