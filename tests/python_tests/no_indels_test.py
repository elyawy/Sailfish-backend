from msasim.simulator import Simulator
from msasim.protocol import SimProtocol
from msasim.constants import MODEL_CODES, ALPHABET_CODES
from msasim.substitutions import ReplacementModelSpec

ROOT_SEQUENCE_LENGTH = 1000

sim_protocol = SimProtocol("tests/trees/normalbranches_nLeaves30000.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)
sim_protocol.set_insertion_rates(0.0)
sim_protocol.set_deletion_rates(0.0)

replacement_model_spec = ReplacementModelSpec(
    alphabet=ALPHABET_CODES.DNA,
    model=MODEL_CODES.NUCJC)

simulation = Simulator(sim_protocol, replacement_model_spec)


msa = simulation()
assert msa.get_length() == ROOT_SEQUENCE_LENGTH

for sparse_seq in msa.get_sparse_msa():
    assert len(sparse_seq) == 1
    assert sparse_seq[0] == ROOT_SEQUENCE_LENGTH

