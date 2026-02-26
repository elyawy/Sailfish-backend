"""Multiple sequence alignment output"""

import _Sailfish

from typing import Dict, List

class Msa:
    """MSA result from simulation"""
    
    def __init__(self, events: list[_Sailfish.IndelEvent], simContext: _Sailfish.SimulationContext):
        self._simContext = simContext
        self._msa = _Sailfish.Msa(events, simContext)

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
    
    def get_msa(self) -> str:
        return self._msa.get_msa_string()
    
    def write_msa(self, file_path) -> None:
        self._msa.write_msa(file_path)

    def get_per_site_rate_categories(self) -> List[int]:
        return self._msa.get_per_site_rate_categories()