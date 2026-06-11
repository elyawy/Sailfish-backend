# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))

from msasim import sailfish as sim

ROOT_SEQUENCE_LENGTH = 100

sim_protocol = sim.SimProtocol("tests/trees/normalbranches_nLeaves10.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)
sim_protocol.set_deletion_rates(0.09)
sim_protocol.set_insertion_rates(0.03)


simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.NOSUBS)
# simulation.set_replacement_model(model=sim.MODEL_CODES.NUCJC)

msa = simulation()

for row in range(msa.get_num_sequences()):
    print(msa.get_msa_row(row))

