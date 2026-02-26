"""Constants and enums for msasim"""

import _Sailfish
from enum import Enum

MODEL_CODES = _Sailfish.modelCode

SITE_RATE_MODELS = _Sailfish.SiteRateModel

class SIMULATION_TYPE(Enum):
    NOSUBS = 0
    DNA = 1
    PROTEIN = 2


SIMULATION_TYPES = set(item for item in SIMULATION_TYPE)