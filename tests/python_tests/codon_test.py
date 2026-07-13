# import sys, pathlib, time
# sys.path.insert(0,str(pathlib.Path(".").resolve()))

from msasim import MODEL_CODES, ReplacementModelSpec, Simulator, SimProtocol, ALPHABET_CODES

ROOT_SEQUENCE_LENGTH = 100

sim_protocol = SimProtocol("tests/trees/normalbranches_nLeaves10.treefile")
sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)
sim_protocol.set_deletion_rates(0.09)
sim_protocol.set_insertion_rates(0.03)

# set a codon model up
replacement_model_spec = ReplacementModelSpec(
    alphabet=ALPHABET_CODES.CODON,
    model=MODEL_CODES.EMPIRICODON)

simulation = Simulator(sim_protocol, replacement_model=replacement_model_spec)

msa = simulation()

for row in range(msa.get_num_sequences()):
    print(msa.get_msa_row(row))

