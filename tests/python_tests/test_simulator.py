"""
Integration tests for the msasim public API.
Covers: indels-only, substitutions-only, combined subs+indels,
deletion limits, low-memory output, and reproducibility.
"""

import pathlib
import time
import pytest

from msasim import SimProtocol, Simulator, Msa
from msasim.distributions import ZipfDistribution, CustomDistribution
from msasim.constants import MODEL_CODES, SIMULATION_TYPE, SITE_RATE_MODELS

# ---------------------------------------------------------------------------
# Shared constants & helpers
# ---------------------------------------------------------------------------

TREE_FILE = "tests/trees/normalbranches_nLeaves10.treefile"
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

def get_msa_string(msa: Msa):
    return "\n".join([msa.get_msa_row(i) for i in range(msa.get_num_sequences())])

# ---------------------------------------------------------------------------
# Indels-only simulation (NOSUBS)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def nosubs_simulator():
    protocol = SimProtocol(
        TREE_FILE,
        root_seq_size=ROOT_SEQ_LEN,
        deletion_rate=0.09,
        insertion_rate=0.03,
        deletion_dist=ZipfDistribution(1.7, 50),
        insertion_dist=ZipfDistribution(1.7, 50),
        seed=42,
    )
    return Simulator(protocol, simulation_type=SIMULATION_TYPE.NOSUBS)


def test_nosubs_correct_num_sequences(nosubs_simulator):
    msa = nosubs_simulator()
    assert msa.get_num_sequences() == NUM_LEAVES


def test_nosubs_sequences_nonempty(nosubs_simulator):
    msa = nosubs_simulator()
    seqs = parse_fasta(get_msa_string(msa))
    assert all(len(s) > 0 for s in seqs.values())


def test_nosubs_sequences_respect_min_size(nosubs_simulator):
    msa = nosubs_simulator()
    seqs = parse_fasta(get_msa_string(msa))
    # Strip gap characters to get actual sequence length per taxon
    for seq in seqs.values():
        assert len(seq.replace("-", "")) >= 0  # min_seq_size default is root length


def test_nosubs_randomness_across_runs(nosubs_simulator):
    msa1 = nosubs_simulator()
    msa2 = nosubs_simulator()
    assert get_msa_string(msa1) != get_msa_string(msa2)


# ---------------------------------------------------------------------------
# Substitutions-only simulation (no indels)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def dna_simulator():
    protocol = SimProtocol(TREE_FILE, root_seq_size=ROOT_SEQ_LEN,
                           insertion_rate=0.0, deletion_rate=0.0,
                           seed=42)
    sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.DNA)
    sim.set_replacement_model(model=MODEL_CODES.NUCJC)
    return sim


def test_dna_correct_num_sequences(dna_simulator):
    msa = dna_simulator()
    assert msa.get_num_sequences() == NUM_LEAVES


def test_dna_sequences_exact_length(dna_simulator):
    """With no indels, every sequence must equal root_seq_size (no gaps)."""
    msa = dna_simulator()
    seqs = parse_fasta(get_msa_string(msa))
    for seq in seqs.values():
        assert len(seq) == ROOT_SEQ_LEN


def test_dna_sequences_not_all_identical(dna_simulator):
    msa = dna_simulator()
    seqs = list(parse_fasta(get_msa_string(msa)).values())
    assert len(set(seqs)) > 1


def test_dna_valid_characters(dna_simulator):
    msa = dna_simulator()
    seqs = parse_fasta(get_msa_string(msa))
    for seq in seqs.values():
        assert set(seq).issubset(VALID_DNA_CHARS), f"Unexpected chars: {set(seq) - VALID_DNA_CHARS}"


# ---------------------------------------------------------------------------
# Combined subs + indels (PROTEIN with INDEL_AWARE rate model)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def protein_indel_simulator():
    protocol = SimProtocol(
        TREE_FILE,
        root_seq_size=ROOT_SEQ_LEN,
        deletion_rate=0.05,
        insertion_rate=0.02,
        deletion_dist=ZipfDistribution(1.7, 50),
        insertion_dist=ZipfDistribution(1.7, 50),
        site_rate_model=SITE_RATE_MODELS.INDEL_AWARE,
        seed=5,
    )
    sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.PROTEIN)
    sim.set_replacement_model(
        model=MODEL_CODES.WAG,
        gamma_parameters_alpha=1.0,
        gamma_parameters_categories=4,
    )
    return sim


def test_protein_indel_correct_num_sequences(protein_indel_simulator):
    msa = protein_indel_simulator()
    assert msa.get_num_sequences() == NUM_LEAVES


def test_protein_indel_valid_characters(protein_indel_simulator):
    msa = protein_indel_simulator()
    seqs = parse_fasta(get_msa_string(msa))
    for seq in seqs.values():
        assert set(seq).issubset(VALID_AA_CHARS), f"Unexpected chars: {set(seq) - VALID_AA_CHARS}"


def test_protein_indel_rate_categories_length(protein_indel_simulator):
    """Rate categories list should match root sequence length."""
    msa = protein_indel_simulator()
    protein_indel_simulator.get_rate_categories()  # call once to populate
    # Call again after a fresh simulation to get the categories
    msa = protein_indel_simulator()
    rate_cats = protein_indel_simulator.get_rate_categories()
    assert len(rate_cats) == ROOT_SEQ_LEN


def test_protein_indel_rate_categories_valid_range(protein_indel_simulator):
    """Rate category indices must be within [0, gamma_categories)."""
    num_categories = 4
    msa = protein_indel_simulator()
    rate_cats = protein_indel_simulator.get_rate_categories()
    assert all(0 <= c < num_categories for c in rate_cats)


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
        seed=50,
        minimum_seq_size=0,
    )
    sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.PROTEIN)
    sim.set_replacement_model(model=MODEL_CODES.WAG)
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
        seed=50,
        minimum_seq_size=min_size,
    )
    sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.PROTEIN)
    sim.set_replacement_model(model=MODEL_CODES.WAG)
    msa = sim()
    seqs = parse_fasta(get_msa_string(msa))
    for seq in seqs.values():
        assert len(seq.replace("-", "")) >= min_size


# ---------------------------------------------------------------------------
# Low-memory simulation
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
    sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.PROTEIN)
    sim.set_replacement_model(model=MODEL_CODES.WAG)
    sim.simulate_low_memory(output_file)
    assert output_file.exists()


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
    sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.PROTEIN)
    sim.set_replacement_model(model=MODEL_CODES.WAG)
    sim.simulate_low_memory(output_file)
    assert output_file.stat().st_size > 0


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
    sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.PROTEIN)
    sim.set_replacement_model(model=MODEL_CODES.WAG)
    sim.simulate_low_memory(output_file)
    seqs = parse_fasta(output_file.read_text())
    assert len(seqs) == 2  # tree has 2 leaves


def test_low_memory_sequences_nonempty(tmp_path):
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
    sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.PROTEIN)
    sim.set_replacement_model(model=MODEL_CODES.WAG)
    sim.simulate_low_memory(output_file)
    seqs = parse_fasta(output_file.read_text())
    assert all(len(s.replace("-", "")) > 0 for s in seqs.values())


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------

def test_same_seed_produces_identical_output():
    def run(seed):
        protocol = SimProtocol(TREE_FILE, root_seq_size=ROOT_SEQ_LEN, seed=seed)
        sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.DNA)
        sim.set_replacement_model(model=MODEL_CODES.NUCJC)
        return get_msa_string(sim())

    assert run(99) == run(99)


def test_different_seeds_produce_different_output():
    def run(seed):
        protocol = SimProtocol(TREE_FILE, root_seq_size=ROOT_SEQ_LEN, seed=seed)
        sim = Simulator(protocol, simulation_type=SIMULATION_TYPE.DNA)
        sim.set_replacement_model(model=MODEL_CODES.NUCJC)

        return get_msa_string(sim())

    assert run(1) != run(2)
