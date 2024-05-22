import sys, pathlib, time
sys.path.insert(0,str(pathlib.Path(".").resolve()))

import numpy as np

from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES

ROOT_SEQUENCE_LENGTH = 1000
NUMBER_OF_MSAS_TO_SIMULATE = 100000

sim_protocol = sim.SimProtocol("tests/trees/normalbranches_nLeaves10.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.DNA)
k2p_params = np.random.uniform(0.1,5, size=NUMBER_OF_MSAS_TO_SIMULATE)

for idx, param in enumerate(k2p_params):
    if idx % 1000 == 0:
        print(f"simulated {idx}\t| {NUMBER_OF_MSAS_TO_SIMULATE}")
    simulation.set_replacement_model(model=MODEL_CODES.HKY,
                                    model_parameters=[0.25,0.25,0.25,0.25, param])
    msa = simulation.simulate()
    # msa.write_msa(f"msas/k2p_param_{param}.fasta")