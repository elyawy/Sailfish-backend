# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))
# Fails currently!
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES

ROOT_SEQUENCE_LENGTH = 10

sim_protocol = sim.SimProtocol("tests/trees/normalbranches_nLeaves10.treefile",
                               deletion_rate=0.1,
                               insertion_rate=0.0,
                               deletion_dist=sim.CustomDistribution([0.0,1.0,0.0]),
                                seed=2)
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.DNA)
# simulation.save_root_sequence()
simulation.set_replacement_model(model=MODEL_CODES.NUCJC, invariant_sites_proportion=0.5)
msa = simulation()

msa.print_msa()