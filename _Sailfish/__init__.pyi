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
import collections.abc
import typing
__all__ = ['AAJC', 'AminoSubstitutionSimulator', 'CPREV45', 'CUSTOM', 'CategorySampler', 'DAYHOFF', 'DELETION', 'DiscreteDistribution', 'EHO_EXTENDED', 'EHO_HELIX', 'EHO_OTHER', 'EMPIRICODON', 'EX_BURIED', 'EX_EHO_BUR_EXT', 'EX_EHO_BUR_HEL', 'EX_EHO_BUR_OTH', 'EX_EHO_EXP_EXT', 'EX_EHO_EXP_HEL', 'EX_EHO_EXP_OTH', 'EX_EXPOSED', 'GTR', 'GammaDistribution', 'HIVB', 'HIVW', 'HKY', 'INDEL_AWARE', 'INSERTION', 'IndelEvent', 'IndelEventType', 'IndelSimulator', 'IntVector', 'JONES', 'LG', 'MTREV24', 'Msa', 'NUCJC', 'NucleotideSubstitutionSimulator', 'SIMPLE', 'SimProtocol', 'SimulationContext', 'SiteRateModel', 'SparseMSA', 'SparseSequenceContainer', 'TAMURA92', 'Tree', 'WAG', 'modelCode', 'modelFactory', 'node']
class AminoSubstitutionSimulator:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: modelFactory, arg1: SimulationContext) -> None:
        ...
    def clear_rates_vec(self) -> None:
        ...
    def get_per_site_rate_categories(self) -> IntVector:
        ...
    def get_site_rates(self) -> list[float]:
        ...
    def init_substitution_sim(self, arg0: modelFactory) -> None:
        ...
    def set_aligned_sequence_map(self, arg0: Msa) -> None:
        ...
    def set_per_site_rate_categories(self, arg0: IntVector) -> None:
        ...
    def set_save_rates(self, arg0: bool) -> None:
        ...
    def simulate_and_write_substitutions(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: str) -> None:
        ...
    def simulate_substitutions(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> SparseSequenceContainer:
        ...
class CategorySampler:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
class DiscreteDistribution:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
class GammaDistribution:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def getAllRates(self) -> list[float]:
        ...
    def getAllRatesProb(self) -> list[float]:
        ...
class IndelEvent:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __repr__(self) -> str:
        ...
    @property
    def length(self) -> int:
        ...
    @property
    def position(self) -> int:
        ...
    @property
    def type(self) -> IndelEventType:
        ...
class IndelEventType:
    """
    Members:
    
      INSERTION
    
      DELETION
    """
    DELETION: typing.ClassVar[IndelEventType]  # value = <IndelEventType.DELETION: 1>
    INSERTION: typing.ClassVar[IndelEventType]  # value = <IndelEventType.INSERTION: 0>
    __members__: typing.ClassVar[dict[str, IndelEventType]]  # value = {'INSERTION': <IndelEventType.INSERTION: 0>, 'DELETION': <IndelEventType.DELETION: 1>}
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
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class IndelSimulator:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: SimulationContext, arg1: SimProtocol) -> None:
        ...
    def generate_events(self) -> list[list[IndelEvent]]:
        ...
    def update_protocol(self, arg0: SimProtocol) -> None:
        ...
class IntVector:
    __hash__: typing.ClassVar[None] = None
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    def __contains__(self, x: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Return true the container contains ``x``
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    def __eq__(self, arg0: IntVector) -> bool:
        ...
    @typing.overload
    def __getitem__(self, s: slice) -> IntVector:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> int:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: IntVector) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[int]:
        ...
    def __len__(self) -> int:
        ...
    def __ne__(self, arg0: IntVector) -> bool:
        ...
    def __repr__(self) -> str:
        """
        Return the canonical string representation of this list.
        """
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: IntVector) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    def count(self, x: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
        Return the number of times ``x`` appears in the list
        """
    @typing.overload
    def extend(self, L: IntVector) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> int:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
        Remove and return the item at index ``i``
        """
    def remove(self, x: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Remove the first item from the list whose value is x. It is an error if there is no such item.
        """
class Msa:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    @typing.overload
    def __init__(self, arg0: collections.abc.Sequence[collections.abc.Sequence[IndelEvent]], arg1: SimulationContext) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: SimulationContext) -> None:
        ...
    def fill_substitutions(self, arg0: SparseSequenceContainer) -> None:
        ...
    def get_msa_row_string(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> str:
        ...
    def get_per_site_rate_categories(self) -> IntVector:
        ...
    def get_sparse_msa(self) -> SparseMSA:
        ...
    def length(self) -> int:
        ...
    def num_sequences(self) -> int:
        ...
    def print_msa(self) -> None:
        ...
    def write_msa(self, arg0: str) -> None:
        ...
class NucleotideSubstitutionSimulator:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: modelFactory, arg1: SimulationContext) -> None:
        ...
    def clear_rates_vec(self) -> None:
        ...
    def get_per_site_rate_categories(self) -> IntVector:
        ...
    def get_site_rates(self) -> list[float]:
        ...
    def init_substitution_sim(self, arg0: modelFactory) -> None:
        ...
    def set_aligned_sequence_map(self, arg0: Msa) -> None:
        ...
    def set_per_site_rate_categories(self, arg0: IntVector) -> None:
        ...
    def set_save_rates(self, arg0: bool) -> None:
        ...
    def simulate_and_write_substitutions(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: str) -> None:
        ...
    def simulate_substitutions(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> SparseSequenceContainer:
        ...
class SimProtocol:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def get_deletion_length_distribution(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> DiscreteDistribution:
        ...
    def get_deletion_rate(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> float:
        ...
    def get_insertion_length_distribution(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> DiscreteDistribution:
        ...
    def get_insertion_rate(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> float:
        ...
    def get_max_insertion_length(self) -> int:
        ...
    def get_minimum_sequence_size(self) -> int:
        ...
    def get_sequence_size(self) -> int:
        ...
    def get_site_rate_model(self) -> SiteRateModel:
        ...
    def set_deletion_length_distributions(self, arg0: collections.abc.Sequence[DiscreteDistribution]) -> None:
        ...
    def set_deletion_rates(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    def set_insertion_length_distributions(self, arg0: collections.abc.Sequence[DiscreteDistribution]) -> None:
        ...
    def set_insertion_rates(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    def set_max_insertion_length(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def set_minimum_sequence_size(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def set_sequence_size(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def set_site_rate_model(self, arg0: SiteRateModel) -> None:
        ...
class SimulationContext:
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __init__(self, tree: Tree, seed: typing.SupportsInt | typing.SupportsIndex, protocol: SimProtocol = None) -> None:
        ...
    def get_category_sampler(self) -> CategorySampler:
        ...
    def get_indel_protocol(self) -> SimProtocol:
        ...
    def get_nodes_to_save(self) -> list[bool]:
        ...
    def get_tree(self) -> Tree:
        ...
    def reseed(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def set_category_sampler(self, arg0: CategorySampler) -> None:
        ...
    def set_protocol(self, arg0: SimProtocol) -> None:
        ...
    def set_save_all(self) -> None:
        ...
    def set_save_leaves(self) -> None:
        ...
    def set_save_root(self) -> None:
        ...
class SiteRateModel:
    """
    Members:
    
      SIMPLE
    
      INDEL_AWARE
    """
    INDEL_AWARE: typing.ClassVar[SiteRateModel]  # value = <SiteRateModel.INDEL_AWARE: 1>
    SIMPLE: typing.ClassVar[SiteRateModel]  # value = <SiteRateModel.SIMPLE: 0>
    __members__: typing.ClassVar[dict[str, SiteRateModel]]  # value = {'SIMPLE': <SiteRateModel.SIMPLE: 0>, 'INDEL_AWARE': <SiteRateModel.INDEL_AWARE: 1>}
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
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class SparseMSA:
    __hash__: typing.ClassVar[None] = None
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    def __contains__(self, x: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> bool:
        """
        Return true the container contains ``x``
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    def __eq__(self, arg0: SparseMSA) -> bool:
        ...
    @typing.overload
    def __getitem__(self, s: slice) -> SparseMSA:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: SparseMSA) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[list[int]]:
        ...
    def __len__(self) -> int:
        ...
    def __ne__(self, arg0: SparseMSA) -> bool:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: SparseMSA) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    def count(self, x: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> int:
        """
        Return the number of times ``x`` appears in the list
        """
    @typing.overload
    def extend(self, L: SparseMSA) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> list[int]:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> list[int]:
        """
        Remove and return the item at index ``i``
        """
    def remove(self, x: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
        Remove the first item from the list whose value is x. It is an error if there is no such item.
        """
class SparseSequenceContainer:
    __hash__: typing.ClassVar[None] = None
    @staticmethod
    def _pybind11_conduit_v1_(*args, **kwargs):
        ...
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    def __contains__(self, x: str) -> bool:
        """
        Return true the container contains ``x``
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    def __eq__(self, arg0: SparseSequenceContainer) -> bool:
        ...
    @typing.overload
    def __getitem__(self, s: slice) -> SparseSequenceContainer:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> str:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: SparseSequenceContainer) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[str]:
        ...
    def __len__(self) -> int:
        ...
    def __ne__(self, arg0: SparseSequenceContainer) -> bool:
        ...
    def __repr__(self) -> str:
        """
        Return the canonical string representation of this list.
        """
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: str) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: SparseSequenceContainer) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: str) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    def count(self, x: str) -> int:
        """
        Return the number of times ``x`` appears in the list
        """
    @typing.overload
    def extend(self, L: SparseSequenceContainer) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: str) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> str:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> str:
        """
        Remove and return the item at index ``i``
        """
    def remove(self, x: str) -> None:
        """
        Remove the first item from the list whose value is x. It is an error if there is no such item.
        """
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
    def root(self) -> node:
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
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __init__(self) -> None:
        ...
    def build_replacement_model(self) -> None:
        ...
    def get_rate_category_sampler(self, max_path_length: typing.SupportsInt | typing.SupportsIndex = 0) -> CategorySampler:
        ...
    def is_model_valid(self) -> bool:
        ...
    def reset(self) -> None:
        ...
    def set_amino_replacement_model_file(self, arg0: str) -> None:
        ...
    def set_model_parameters(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    def set_replacement_model(self, arg0: modelCode) -> None:
        ...
    def set_site_rate_model(self, rates: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], stationary_probs: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], transition_matrix: collections.abc.Sequence[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]] = []) -> None:
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
AAJC: modelCode  # value = <modelCode.AAJC: 1>
CPREV45: modelCode  # value = <modelCode.CPREV45: 6>
CUSTOM: modelCode  # value = <modelCode.CUSTOM: 26>
DAYHOFF: modelCode  # value = <modelCode.DAYHOFF: 7>
DELETION: IndelEventType  # value = <IndelEventType.DELETION: 1>
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
INDEL_AWARE: SiteRateModel  # value = <SiteRateModel.INDEL_AWARE: 1>
INSERTION: IndelEventType  # value = <IndelEventType.INSERTION: 0>
JONES: modelCode  # value = <modelCode.JONES: 8>
LG: modelCode  # value = <modelCode.LG: 13>
MTREV24: modelCode  # value = <modelCode.MTREV24: 9>
NUCJC: modelCode  # value = <modelCode.NUCJC: 0>
SIMPLE: SiteRateModel  # value = <SiteRateModel.SIMPLE: 0>
TAMURA92: modelCode  # value = <modelCode.TAMURA92: 4>
WAG: modelCode  # value = <modelCode.WAG: 10>
