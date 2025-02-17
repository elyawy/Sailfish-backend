# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))

from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES, CustomDistribution

ROOT_SEQUENCE_LENGTH = 100


sim_protocol = sim.SimProtocol("(A:0.5,B:0.5);",
                               deletion_rate=11.0,
                               insertion_rate=0.00,
                               deletion_dist=CustomDistribution([1.0]),
                               insertion_dist=CustomDistribution([1.0]),
                               seed=50, minimum_seq_size=0)
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)

simulation.set_replacement_model(model=MODEL_CODES.WAG, 
                                 gamma_parameters_alpha=1.0, 
                                 gamma_parameters_categories=4)
# simulation.save_root_sequence()
msa = simulation()
msa.print_msa()
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