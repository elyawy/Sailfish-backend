# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))

from msasim import Simulator, SimProtocol

ROOT_SEQUENCE_LENGTH = 100

sim_protocol = SimProtocol("tests/trees/normalbranches_nLeaves10.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)
sim_protocol.set_deletion_rates(0.09)
sim_protocol.set_insertion_rates(0.03)


simulation = Simulator(sim_protocol)
# simulation.set_replacement_model(model=sim.MODEL_CODES.NUCJC)

msa = simulation()

for row in range(msa.get_num_sequences()):
    print(msa.get_msa_row(row))

