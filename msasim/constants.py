"""Constants and enums for msasim"""

import _Sailfish
from enum import Enum

MODEL_CODES = _Sailfish.modelCode

class SIMULATION_TYPE(Enum):
    NOSUBS = 0
    DNA = 1
    PROTEIN = 2