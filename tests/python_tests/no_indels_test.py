# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))

from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES

ROOT_SEQUENCE_LENGTH = 1000

sim_protocol = sim.SimProtocol("tests/trees/normalbranches_nLeaves30000.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.DNA)

simulation.set_replacement_model(model=MODEL_CODES.NUCJC)
msa = simulation()

msa.print_msa()