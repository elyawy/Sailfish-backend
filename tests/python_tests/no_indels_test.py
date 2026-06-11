from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES

ROOT_SEQUENCE_LENGTH = 1000

sim_protocol = sim.SimProtocol("tests/trees/normalbranches_nLeaves30000.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)
sim_protocol.set_insertion_rates(0.0)
sim_protocol.set_deletion_rates(0.0)
simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.DNA)

simulation.set_replacement_model(model=MODEL_CODES.NUCJC)
msa = simulation()
assert msa.get_length() == ROOT_SEQUENCE_LENGTH

for sparse_seq in msa.get_sparse_msa():
    assert len(sparse_seq) == 1
    assert sparse_seq[0] == ROOT_SEQUENCE_LENGTH

