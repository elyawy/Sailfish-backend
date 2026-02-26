# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))

from msasim import SimProtocol
from msasim import Simulator
from msasim.distributions import ZipfDistribution
from msasim.constants import MODEL_CODES, SIMULATION_TYPE, SITE_RATE_MODELS


ROOT_SEQUENCE_LENGTH = 100

sim_protocol = SimProtocol("tests/trees/normalbranches_nLeaves10.treefile",
                               deletion_rate=0.0,
                               insertion_rate=0.02,
                               deletion_dist=ZipfDistribution(1.7, 50),
                               insertion_dist=ZipfDistribution(1.7, 50),
                               site_rate_model=SITE_RATE_MODELS.INDEL_AWARE,
                               seed=5)
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

simulation = Simulator(sim_protocol, simulation_type=SIMULATION_TYPE.PROTEIN)
# simulation.save_root_sequence()
simulation.set_replacement_model(model=MODEL_CODES.WAG, 
                                 gamma_parameters_alpha=1.0, 
                                 gamma_parameters_categories=4)


msa = simulation()
msa.print_msa()

rate_categories =simulation.get_rate_categories()

print("".join(map(str, rate_categories)))
# counter = 0
# while True:
#     print(counter)
#     msa = simulation()
#     counter += 1

# for i in range(10000):
#     print(i)
#     msa = simulation()
# msa.print_msa()
# msa.print_msa()