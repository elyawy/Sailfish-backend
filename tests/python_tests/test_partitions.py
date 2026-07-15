"""Basic tests for the Partitions / PartitionModel composite API."""

import pytest

from msasim import  SimProtocol
from msasim.substitutions import ReplacementModelSpec
from msasim.advanced import PartitionModel, Partitions
from msasim.constants import MODEL_CODES, ALPHABET_CODES


def parse_fasta(fasta_str: str) -> dict:
    sequences = {}
    current_name = None
    for line in fasta_str.strip().splitlines():
        if line.startswith(">"):
            current_name = line[1:].strip()
            sequences[current_name] = ""
        elif current_name is not None:
            sequences[current_name] += line.strip()
    return sequences


def test_partitions_same_taxa_concatenates_lengths(tmp_path):
    """Two partitions on the same tree -> every taxon's sequence length is the sum of both."""
    output_file = tmp_path / "partitions.fasta"

    partition1 = PartitionModel(
        name="p1",
        replacement_model=ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA),
        indel_model=SimProtocol(
            "(A:0.1,B:0.1);", root_seq_size=50,
            insertion_rate=0.0, deletion_rate=0.0, seed=1,
        ),
    )
    partition2 = PartitionModel(
        name="p2",
        replacement_model=ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA),
        indel_model=SimProtocol(
            "(A:0.1,B:0.1);", root_seq_size=30,
            insertion_rate=0.0, deletion_rate=0.0, seed=2,
        ),
    )

    partitions = Partitions()
    partitions.add_partition(partition1)
    partitions.add_partition(partition2)
    partitions.simulate(output_file)

    seqs = parse_fasta(output_file.read_text())
    assert set(seqs.keys()) == {"A", "B"}
    for seq in seqs.values():
        assert len(seq) == 50 + 30


def test_partitions_different_taxa_pads_with_gaps(tmp_path):
    """Taxon missing from one partition is padded with gaps of that partition's length."""
    output_file = tmp_path / "partitions.fasta"

    partition1 = PartitionModel(
        name="p1",
        replacement_model=ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA),
        indel_model=SimProtocol(
            "(A:0.1,B:0.1);", root_seq_size=20,
            insertion_rate=0.0, deletion_rate=0.0, seed=1,
        ),
    )
    partition2 = PartitionModel(
        name="p2",
        replacement_model=ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA),
        indel_model=SimProtocol(
            "(A:0.1,C:0.1);", root_seq_size=15,
            insertion_rate=0.0, deletion_rate=0.0, seed=2,
        ),
    )

    partitions = Partitions()
    partitions.add_partition(partition1)
    partitions.add_partition(partition2)
    partitions.simulate(output_file)

    seqs = parse_fasta(output_file.read_text())
    assert set(seqs.keys()) == {"A", "B", "C"}
    assert len(seqs["A"]) == 35
    assert len(seqs["B"]) == 35
    assert seqs["B"][20:] == "-" * 15  # B missing from partition 2
    assert len(seqs["C"]) == 35
    assert seqs["C"][:20] == "-" * 20  # C missing from partition 1


def test_partitions_no_partitions_raises():
    partitions = Partitions()
    with pytest.raises(ValueError):
        partitions.simulate("output.fasta")


def test_partitions_visual_inspection(tmp_path):
    """Not a real assertion test -- run with `pytest -s` to print and eyeball the merged FASTA."""
    output_file = tmp_path / "partitions_visual.fasta"

    partition1 = PartitionModel(
        name="p1",
        replacement_model=ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA),
        indel_model=SimProtocol(
            "(A:0.1,B:0.1);", root_seq_size=15,
            insertion_rate=0.0, deletion_rate=0.0, seed=1,
        ),
    )
    partition2 = PartitionModel(
        name="p2",
        replacement_model=ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA),
        indel_model=SimProtocol(
            "(A:0.1,C:0.1);", root_seq_size=10,
            insertion_rate=0.0, deletion_rate=0.0, seed=2,
        ),
    )

    partitions = Partitions()
    partitions.add_partition(partition1)
    partitions.add_partition(partition2)
    partitions.simulate(output_file)

    content = output_file.read_text()
    print("\n--- merged partitions FASTA ---")
    print(content)
    print("--- end ---")