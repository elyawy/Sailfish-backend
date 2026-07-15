"""
Integration tests for the msasim public API.
Covers: indels-only, subs-only, full (indels+subs), rate tracking,
sparse MSA, deletion limits, low-memory output, and reproducibility.
"""

import pathlib
import pytest
import warnings

from msasim import SimProtocol, Simulator, ReplacementModelSpec, SiteRateModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES
from msasim.distributions import ZipfDistribution, CustomDistribution

# ---------------------------------------------------------------------------
# Shared constants & helpers
# ---------------------------------------------------------------------------

TREE_FILE = str(pathlib.Path(__file__).parent.parent / "trees" / "normalbranches_nLeaves10.treefile")
NUM_LEAVES = 10
ROOT_SEQ_LEN = 100

VALID_DNA_CHARS = set("ACGT-")
VALID_AA_CHARS = set("ACDEFGHIKLMNPQRSTVWY-")


def parse_fasta(fasta_str: str) -> dict:
    """Return {seq_name: sequence_string} from a FASTA-format string."""
    sequences = {}
    current_name = None
    for line in fasta_str.strip().splitlines():
        if line.startswith(">"):
            current_name = line[1:].strip()
            sequences[current_name] = ""
        elif current_name is not None:
            sequences[current_name] += line.strip()
    return sequences


def get_msa_rows(msa) -> list:
    return [msa.get_msa_row(i) for i in range(msa.get_num_sequences())]


def get_sequences(msa) -> dict:
    """Parse all rows into {name: sequence} dict."""
    result = {}
    for row in get_msa_rows(msa):
        name, seq = row.split("\n", 1)
        result[name.lstrip(">")] = seq
    return result


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def nosubs_simulator():
    """Indels only — no replacement_model passed."""
    protocol = SimProtocol(
        TREE_FILE,
        root_seq_size=ROOT_SEQ_LEN,
        deletion_rate=0.09,
        insertion_rate=0.03,
        deletion_dist=ZipfDistribution(1.7, 50),
        insertion_dist=ZipfDistribution(1.7, 50),
        seed=42,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Simulator(protocol)


@pytest.fixture(scope="module")
def dna_simulator():
    """DNA substitutions only — zero indel rates."""
    protocol = SimProtocol(
        TREE_FILE,
        root_seq_size=ROOT_SEQ_LEN,
        insertion_rate=0.0,
        deletion_rate=0.0,
        seed=42,
    )
    rep_model = ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA)
    return Simulator(protocol, replacement_model=rep_model)


@pytest.fixture(scope="module")
def protein_indel_simulator():
    """Protein + indels with gamma rate heterogeneity."""
    rate_model = SiteRateModelSpec(
        gamma_alpha=1.0,
        gamma_categories=4,
    )
    rep_model = ReplacementModelSpec(
        model=MODEL_CODES.WAG,
        alphabet=ALPHABET_CODES.PROTEIN,
        site_rate_model=rate_model,
    )
    protocol = SimProtocol(
        TREE_FILE,
        root_seq_size=ROOT_SEQ_LEN,
        deletion_rate=0.05,
        insertion_rate=0.02,
        deletion_dist=ZipfDistribution(1.7, 50),
        insertion_dist=ZipfDistribution(1.7, 50),
        seed=5,
    )
    return Simulator(protocol, replacement_model=rep_model)


# ---------------------------------------------------------------------------
# Indels-only (NOSUBS)
# ---------------------------------------------------------------------------

def test_nosubs_correct_num_sequences(nosubs_simulator):
    msa = nosubs_simulator()
    assert msa.get_num_sequences() == NUM_LEAVES


def test_nosubs_sequences_nonempty(nosubs_simulator):
    msa = nosubs_simulator()
    seqs = get_sequences(msa)
    assert all(len(s) > 0 for s in seqs.values())


def test_nosubs_randomness_across_runs(nosubs_simulator):
    msa1 = nosubs_simulator()
    msa2 = nosubs_simulator()
    assert get_msa_rows(msa1) != get_msa_rows(msa2)


# ---------------------------------------------------------------------------
# Substitutions-only (DNA, zero indel rates)
# ---------------------------------------------------------------------------

def test_dna_correct_num_sequences(dna_simulator):
    msa = dna_simulator()
    assert msa.get_num_sequences() == NUM_LEAVES


def test_dna_sequences_exact_length(dna_simulator):
    """With no indels, alignment length must equal root_seq_size."""
    msa = dna_simulator()
    assert msa.get_length() == ROOT_SEQ_LEN


def test_dna_sequences_not_all_identical(dna_simulator):
    msa = dna_simulator()
    seqs = list(get_sequences(msa).values())
    assert len(set(seqs)) > 1


def test_dna_valid_characters(dna_simulator):
    msa = dna_simulator()
    for seq in get_sequences(msa).values():
        assert set(seq).issubset(VALID_DNA_CHARS), f"Unexpected chars: {set(seq) - VALID_DNA_CHARS}"


def test_dna_sparse_msa_no_indels(dna_simulator):
    """With no indels, every sparse block should be a single run of ROOT_SEQ_LEN."""
    msa = dna_simulator()
    for sparse_seq in msa.get_sparse_msa():
        assert len(sparse_seq) == 1
        assert sparse_seq[0] == ROOT_SEQ_LEN


# ---------------------------------------------------------------------------
# Full simulation: protein + indels + INDEL_AWARE rate model
# ---------------------------------------------------------------------------

def test_protein_indel_correct_num_sequences(protein_indel_simulator):
    msa = protein_indel_simulator()
    assert msa.get_num_sequences() == NUM_LEAVES


def test_protein_indel_valid_characters(protein_indel_simulator):
    msa = protein_indel_simulator()
    for seq in get_sequences(msa).values():
        assert set(seq).issubset(VALID_AA_CHARS), f"Unexpected chars: {set(seq) - VALID_AA_CHARS}"


def test_protein_indel_rate_categories_length(protein_indel_simulator):
    """Rate categories list should match alignment length (includes inserted sites)."""
    msa = protein_indel_simulator()
    rate_cats = protein_indel_simulator.get_rate_categories()
    assert len(rate_cats) == msa.get_length()


def test_protein_indel_rate_categories_valid_range(protein_indel_simulator):
    """Rate category indices must be within [0, gamma_categories)."""
    msa = protein_indel_simulator()
    rate_cats = protein_indel_simulator.get_rate_categories()
    assert all(0 <= c < 4 for c in rate_cats)


# ---------------------------------------------------------------------------
# Per-site rate tracking (save_rates / get_rates)
# ---------------------------------------------------------------------------

def test_save_rates_returns_correct_length():
    protocol = SimProtocol(
        TREE_FILE,
        root_seq_size=ROOT_SEQ_LEN,
        insertion_rate=0.0,
        deletion_rate=0.0,
        seed=7,
    )
    rate_model = SiteRateModelSpec(gamma_alpha=0.5, gamma_categories=8)
    rep_model = ReplacementModelSpec(
        model=MODEL_CODES.NUCJC,
        alphabet=ALPHABET_CODES.DNA,
        site_rate_model=rate_model,
    )
    sim = Simulator(protocol, replacement_model=rep_model)
    sim.save_rates(True)
    sim()
    rates = sim.get_rates()
    assert len(rates) == ROOT_SEQ_LEN


def test_save_rates_returns_positive_floats():
    protocol = SimProtocol(
        TREE_FILE,
        root_seq_size=ROOT_SEQ_LEN,
        insertion_rate=0.0,
        deletion_rate=0.0,
        seed=8,
    )
    rate_model = SiteRateModelSpec(gamma_alpha=0.5, gamma_categories=4)
    rep_model = ReplacementModelSpec(
        model=MODEL_CODES.NUCJC,
        alphabet=ALPHABET_CODES.DNA,
        site_rate_model=rate_model,
    )
    sim = Simulator(protocol, replacement_model=rep_model)
    sim.save_rates(True)
    sim()
    rates = sim.get_rates()
    assert all(isinstance(r, float) and r > 0 for r in rates)


def test_site_rate_correlation_runs():
    """site_rate_correlation > 0 should not crash."""
    protocol = SimProtocol(
        TREE_FILE,
        root_seq_size=ROOT_SEQ_LEN,
        insertion_rate=0.0,
        deletion_rate=0.0,
        seed=9,
    )
    rate_model = SiteRateModelSpec(
        gamma_alpha=0.5,
        gamma_categories=8,
        site_rate_correlation=0.99,
    )
    rep_model = ReplacementModelSpec(
        model=MODEL_CODES.NUCJC,
        alphabet=ALPHABET_CODES.DNA,
        site_rate_model=rate_model,
    )
    sim = Simulator(protocol, replacement_model=rep_model)
    msa = sim()
    assert msa.get_num_sequences() == NUM_LEAVES




# ---------------------------------------------------------------------------
# Deletion limit / minimum sequence size
# ---------------------------------------------------------------------------

def test_deletion_limit_sequences_exist():
    """Even under extreme deletion pressure, simulation should complete."""
    protocol = SimProtocol(
        "(A:0.5,B:0.5);",
        root_seq_size=ROOT_SEQ_LEN,
        deletion_rate=11.0,
        insertion_rate=0.0,
        deletion_dist=CustomDistribution([1.0]),
        insertion_dist=CustomDistribution([1.0]),
        minimum_seq_size=0,
        seed=50,
    )
    rep_model = ReplacementModelSpec(model=MODEL_CODES.WAG, alphabet=ALPHABET_CODES.PROTEIN)
    sim = Simulator(protocol, replacement_model=rep_model)
    msa = sim()
    assert msa.get_num_sequences() == 2


def test_minimum_seq_size_respected():
    """No sequence should fall below minimum_seq_size under high deletion pressure."""
    min_size = 10
    protocol = SimProtocol(
        "(A:0.5,B:0.5);",
        root_seq_size=ROOT_SEQ_LEN,
        deletion_rate=11.0,
        insertion_rate=0.0,
        deletion_dist=CustomDistribution([1.0]),
        insertion_dist=CustomDistribution([1.0]),
        minimum_seq_size=min_size,
        seed=50,
    )
    rep_model = ReplacementModelSpec(model=MODEL_CODES.WAG, alphabet=ALPHABET_CODES.PROTEIN)
    sim = Simulator(protocol, replacement_model=rep_model)
    msa = sim()
    for seq in get_sequences(msa).values():
        assert len(seq.replace("-", "")) >= min_size


# ---------------------------------------------------------------------------
# Low-memory / output_path mode
# ---------------------------------------------------------------------------

def test_low_memory_file_created(tmp_path):
    output_file = tmp_path / "output.fasta"
    protocol = SimProtocol(
        "(A:0.5,B:0.5);",
        root_seq_size=ROOT_SEQ_LEN,
        deletion_rate=0.01,
        insertion_rate=0.01,
        deletion_dist=ZipfDistribution(1.08, 50),
        insertion_dist=ZipfDistribution(1.08, 50),
        seed=1234,
    )
    rep_model = ReplacementModelSpec(model=MODEL_CODES.WAG, alphabet=ALPHABET_CODES.PROTEIN)
    sim = Simulator(protocol, replacement_model=rep_model)
    sim.simulate(output_path=output_file)
    assert (tmp_path / "output_replicate_1.fasta").exists()


def test_low_memory_file_nonempty(tmp_path):
    output_file = tmp_path / "output.fasta"
    protocol = SimProtocol(
        "(A:0.5,B:0.5);",
        root_seq_size=ROOT_SEQ_LEN,
        deletion_rate=0.01,
        insertion_rate=0.01,
        deletion_dist=ZipfDistribution(1.08, 50),
        insertion_dist=ZipfDistribution(1.08, 50),
        seed=1234,
    )
    rep_model = ReplacementModelSpec(model=MODEL_CODES.WAG, alphabet=ALPHABET_CODES.PROTEIN)
    sim = Simulator(protocol, replacement_model=rep_model)
    sim.simulate(output_path=output_file)
    out = tmp_path / "output_replicate_1.fasta"
    assert out.stat().st_size > 0


def test_low_memory_correct_num_sequences(tmp_path):
    output_file = tmp_path / "output.fasta"
    protocol = SimProtocol(
        "(A:0.5,B:0.5);",
        root_seq_size=ROOT_SEQ_LEN,
        deletion_rate=0.01,
        insertion_rate=0.01,
        deletion_dist=ZipfDistribution(1.08, 50),
        insertion_dist=ZipfDistribution(1.08, 50),
        seed=1234,
    )
    rep_model = ReplacementModelSpec(model=MODEL_CODES.WAG, alphabet=ALPHABET_CODES.PROTEIN)
    sim = Simulator(protocol, replacement_model=rep_model)
    sim.simulate(output_path=output_file)
    out = tmp_path / "output_replicate_1.fasta"
    seqs = parse_fasta(out.read_text())
    assert len(seqs) == 2  # tree has 2 leaves


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

def test_same_seed_produces_identical_output():
    def run(seed):
        protocol = SimProtocol(TREE_FILE, root_seq_size=ROOT_SEQ_LEN, seed=seed)
        rep_model = ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA)
        sim = Simulator(protocol, replacement_model=rep_model)
        return get_msa_rows(sim())

    assert run(99) == run(99)


def test_different_seeds_produce_different_output():
    def run(seed):
        protocol = SimProtocol(TREE_FILE, root_seq_size=ROOT_SEQ_LEN, seed=seed)
        rep_model = ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA)
        sim = Simulator(protocol, replacement_model=rep_model)
        return get_msa_rows(sim())

    assert run(1) != run(2)


# ---------------------------------------------------------------------------
# Codon
# ---------------------------------------------------------------------------

def test_codon_basic():
    protocol = SimProtocol(
        TREE_FILE,
        root_seq_size=ROOT_SEQ_LEN,
        insertion_rate=0.0,
        deletion_rate=0.0,
        seed=42,
    )
    rep_model = ReplacementModelSpec(model=MODEL_CODES.EMPIRICODON, alphabet=ALPHABET_CODES.CODON)
    sim = Simulator(protocol, replacement_model=rep_model)
    msa = sim()
    assert msa.get_num_sequences() == NUM_LEAVES
    for seq in get_sequences(msa).values():
        assert len(seq) == ROOT_SEQ_LEN * 3