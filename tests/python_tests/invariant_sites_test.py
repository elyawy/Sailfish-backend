# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))
# Fails currently!
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES

ROOT_SEQUENCE_LENGTH = 50

sim_protocol = sim.SimProtocol("tests/trees/normalbranches_nLeaves10.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.DNA)

simulation.set_replacement_model(model=MODEL_CODES.NUCJC,gamma_parameters_alpha=1.0,
                                 gamma_parameters_catergories=4, invariant_sites_proportion=0.2)
msa = simulation()

msa.print_msa()