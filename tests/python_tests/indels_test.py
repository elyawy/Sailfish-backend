# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))

from msasim import SimProtocol
from msasim import Simulator
from msasim.distributions import ZipfDistribution
from msasim.constants import MODEL_CODES, SITE_RATE_MODELS
from msasim.substitutions import ReplacementModelSpec, SiteRateModelSpec

ROOT_SEQUENCE_LENGTH = 100

sim_protocol = SimProtocol("tests/trees/normalbranches_nLeaves10.treefile",
                               deletion_rate=0.02,
                               insertion_rate=0.02,
                               deletion_dist=ZipfDistribution(1.7, 50),
                               insertion_dist=ZipfDistribution(1.7, 50),
                               site_rate_model=SITE_RATE_MODELS.SIMPLE,
                               seed=11)
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

rate_model = SiteRateModelSpec(gamma_alpha=1.0, gamma_categories=4)

replacement_model = ReplacementModelSpec(model=MODEL_CODES.WAG, 
                                         site_rate_model=rate_model)

simulator = Simulator(sim_protocol, replacement_model)

simulator.save_root_sequence()
simulator.set_root_sequence("".join(["M" for _ in range(ROOT_SEQUENCE_LENGTH)]))

msa = simulator()
msa.print_msa()
