import numpy as np
from msasim import SimProtocol, MODEL_CODES, ALPHABET_CODES
from msasim.substitutions import ReplacementModelSpec
from msasim.advanced import MixtureModel, Mixture

tree = "(A:0.1,B:0.1,C:0.1);"
models = [ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA) for _ in range(3)]
weights = [0.3, 0.3, 0.4]


def make_mixture(seed, insertion_rate=0.0, deletion_rate=0.0):
    indel_protocol = SimProtocol(
        tree, root_seq_size=50,
        insertion_rate=insertion_rate, deletion_rate=deletion_rate, seed=seed,
    )
    mixture_model = MixtureModel(
        name="test_mixture",
        replacement_models=models,
        model_weights=weights,
        indel_model=indel_protocol,
    )
    return Mixture(mixture_model)


def test_mixture_basic_smoke():
    msa = make_mixture(seed=1).simulate()
    assert msa.get_num_sequences() == 3


def test_mixture_with_indels_runs_and_matches_length():
    """Nonzero indel rate -> exercise gap-stripping path, check no crash."""
    msa = make_mixture(seed=2, insertion_rate=0.05, deletion_rate=0.05).simulate()

    assert msa.get_num_sequences() == 3
    for i in range(msa.get_num_sequences()):
        row = msa.get_msa_row(i).split("\n", 1)[1]
        assert len(row) == msa.get_length()


def test_mixture_same_seed_reproducible():
    """Same seed, fresh objects -> identical output across runs."""
    msa1 = make_mixture(seed=3).simulate()
    msa2 = make_mixture(seed=3).simulate()

    rows1 = [msa1.get_msa_row(i) for i in range(msa1.get_num_sequences())]
    rows2 = [msa2.get_msa_row(i) for i in range(msa2.get_num_sequences())]
    assert rows1 == rows2