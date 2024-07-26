# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))

from msasim import sailfish as sim

ROOT_SEQUENCE_LENGTH = 100

sim_protocol = sim.SimProtocol("tests/trees/normalbranches_nLeaves10.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)
sim_protocol.set_deletion_rates(0.09)
sim_protocol.set_insertion_rates(0.03)

simulation = sim.Simulator(sim_protocol)

for i in range(10):
    msa = simulation()

print(msa.get_msa())