import sys, pathlib, time
sys.path.insert(0,str(pathlib.Path(".").resolve()))

from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES, ZipfDistribution

ROOT_SEQUENCE_LENGTH = 150

sim_protocol = sim.SimProtocol("((A:0.1,B:1.0):0.05,(C:0.2):0.01);",
                               deletion_rate=0.01,
                               insertion_rate=0.01,
                               deletion_dist=ZipfDistribution(1.08, 50),
                               insertion_dist=ZipfDistribution(1.08, 50),
                               seed=10)
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)

simulation.set_replacement_model(model=MODEL_CODES.WAG, 
                                 gamma_parameters_alpha=1.0, 
                                 gamma_parameters_categories=4)
simulation.save_root_sequence()

msa = simulation()
msa.print_msa()
