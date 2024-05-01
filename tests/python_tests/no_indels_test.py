from Sailfish import simulator as sim
from Sailfish.simulator import MODEL_CODES

ROOT_SEQUENCE_LENGTH = 1000

sim_protocol = sim.SimProtocol("tests/trees/normalbranches_nLeaves30000.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.DNA)

simulation.set_replacement_model(model=MODEL_CODES.NUCJC)
msa = simulation()

# msa.print_msa()