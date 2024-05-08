import sys, pathlib, time
sys.path.insert(0,str(pathlib.Path(".").resolve()))

from Sailfish import simulator as sim
from Sailfish.simulator import MODEL_CODES

ROOT_SEQUENCE_LENGTH = 100

sim_protocol = sim.SimProtocol("tests/trees/normalbranches_nLeaves10.treefile",
                               deletion_rate=0.01,
                               insertion_rate=0.01)
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)

simulation.set_replacement_model(model=MODEL_CODES.JONES, gamma_parameters_alpha=1.0, gamma_parameters_catergories=4)
msa = simulation()
msa.print_msa()
# msa.print_msa()