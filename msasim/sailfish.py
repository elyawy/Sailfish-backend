"""
DEPRECATED: Import from msasim directly instead.

This module exists for backwards compatibility only.
Instead of:
    from msasim import sailfish as sim
    sim.Simulator(...)

Use:
    from msasim import Simulator
    Simulator(...)
"""

import warnings
warnings.warn(
    "Importing from msasim.sailfish is deprecated. "
    "Import directly from msasim instead: 'from msasim import Simulator'",
    DeprecationWarning,
    stacklevel=2
)

# Re-export everything for backwards compatibility
from .distributions import *
from .tree import *
from .protocol import *
from .simulator import *
from .msa import *
from .constants import *
