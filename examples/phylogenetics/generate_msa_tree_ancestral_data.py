"""
Generate training data for one of three tasks (chosen via --data_type):

  MSA                 - reconstruct the full (column-major) multiple sequence alignment
  ANCESTRAL_SEQUENCE  - reconstruct the root (ancestral) sequence
  TREE                - reconstruct the guide tree in newick format

For each simulation, a random (or fixed) tree, root length, and indel rates/parameters
are sampled, Sailfish generates the alignment, and one output line is written as:

    <source unaligned sequences><special_token><target>

A special-tokens file is also written, listing every leaf name, common newick
characters, and the three special task tokens, for downstream tokenizer setup.
"""

import argparse
import os
import random

import numpy as np
from ete3 import Tree
from msasim import (
    SimProtocol,
    Simulator,
    ReplacementModelSpec,
    SiteRateModelSpec,
    MODEL_CODES,
    ALPHABET_CODES,
    ZipfDistribution,
)

# task types (--data_type)
MSA = "MSA"
TREE = "TREE"
ANCESTRAL_SEQUENCE = "ANCESTRAL_SEQUENCE"

# special tokens marking the start of the target segment, one per task
SPECIAL_MSA = "<ALIGN>"
SPECIAL_TREE = "<INFER>"
SPECIAL_ANCESTRAL = "<RECONSTRUCT>"

# ete3's default label for the internal node saved as the "root" ancestral sequence
ANCESTRAL_NODE_NAME = "N1"

# hardcoded topology used when --fixed_tree is set (only defined for 8 species)
FIXED_TREE_8_SPECIES = (
    "((l5:0.24,(l2:0.361,l4:0.83):0.954):0.804,"
    "((l1:0.973,(l7:0.016,l0:0.418):0.501):0.083,(l6:0.331,l3:0.601):0.882):0.94);"
)

SPECIAL_CHARACTERS = [
    "(", ")", ",", ".", "/", "\\", "'", '"',
    "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
    ":", "#", "$", "^", "*", "!", "_", ";",
]


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate multi-task simulated training data (MSA / ancestral sequence / tree)"
    )

    parser.add_argument(
        "--data_type", type=str, choices=[ANCESTRAL_SEQUENCE, MSA, TREE], default=MSA,
        help="which target task to generate data for",
    )

    # files
    parser.add_argument("--output_path_folder", type=str, required=True, help="output folder")
    parser.add_argument("--output_file_name", type=str, required=True, help="output file name (source+target pairs)")
    parser.add_argument("--output_source_tokenize_name", type=str, default=None,
                        help="optional file with source sequences only, for tokenizer training")
    parser.add_argument("--output_target_tokenize_name", type=str, default=None,
                        help="optional file with target sequences only, for tokenizer training")
    parser.add_argument("--special_tokens_file_name", type=str, default="special_tokens.txt",
                        help="output file listing special tokens")

    # tree
    parser.add_argument("--num_species", type=int, default=10, help="number of species (leaves)")
    parser.add_argument("--min_branch_length", type=float, default=0.00, help="minimum branch length")
    parser.add_argument("--max_branch_length", type=float, default=1.00, help="maximum branch length")
    parser.add_argument("--fixed_tree", default=False, action=argparse.BooleanOptionalAction,
                        help="use a fixed topology instead of a random one (only defined for num_species=8)")

    # simulator settings
    parser.add_argument("--num_of_simulations", type=int, required=True, help="number of simulations to run")
    parser.add_argument("--min_deletion_rate", type=float, default=0.00, help="minimum deletion rate")
    parser.add_argument("--max_deletion_rate", type=float, default=0.05, help="maximum deletion rate")
    parser.add_argument("--min_insertion_rate", type=float, default=0.00, help="minimum insertion rate")
    parser.add_argument("--max_insertion_rate", type=float, default=0.05, help="maximum insertion rate")
    parser.add_argument("--average_root_len", type=int, default=345, help="mean of root sequence length")
    parser.add_argument("--std_root_len", type=int, default=76, help="std of root sequence length")
    parser.add_argument("--min_root_len", type=int, default=25, help="lower bound on root sequence length")
    parser.add_argument("--min_deletion_a_val", type=float, default=1.00, help="minimum Zipf 'a' parameter for deletions")
    parser.add_argument("--max_deletion_a_val", type=float, default=2.00, help="maximum Zipf 'a' parameter for deletions")
    parser.add_argument("--min_insertion_a_val", type=float, default=1.00, help="minimum Zipf 'a' parameter for insertions")
    parser.add_argument("--max_insertion_a_val", type=float, default=2.00, help="maximum Zipf 'a' parameter for insertions")

    # output formatting
    parser.add_argument("--source_delimiter", type=str, default="|", help="delimiter between sequences in the source")
    parser.add_argument("--write_every", type=int, default=1000, help="flush output every N simulations")

    return parser.parse_args()


class BufferedWriter:
    """Accumulates lines and flushes them to a file every `write_every` lines. No-op if path is None."""

    def __init__(self, path, write_every):
        self.path = path
        self.write_every = write_every
        self.buffer = []

    def add(self, line):
        if self.path is None:
            return
        self.buffer.append(line)
        if len(self.buffer) >= self.write_every:
            self.flush()

    def flush(self):
        if not self.path or not self.buffer:
            return
        with open(self.path, "a") as f:
            for line in self.buffer:
                f.write(line + "\n")
        self.buffer = []


def build_special_tokens_list(leaf_names):
    return leaf_names + SPECIAL_CHARACTERS + [SPECIAL_MSA, SPECIAL_ANCESTRAL, SPECIAL_TREE]


def build_tree(args, leaf_names):
    """Return a newick tree string: either the fixed 8-species topology, or a random one."""
    if args.fixed_tree:
        if args.num_species != 8:
            exit("fixed tree is not defined for this number of species")
        return FIXED_TREE_8_SPECIES

    t = Tree()
    shuffled_names = leaf_names.copy()
    random.shuffle(shuffled_names)
    t.populate(len(leaf_names), names_library=shuffled_names)

    newick = t.write(format=5).replace(":1", "##")
    while "##" in newick:
        branch_length = round(random.uniform(args.min_branch_length, args.max_branch_length), 3)
        newick = newick.replace("##", f":{branch_length}", 1)
    return newick


def sample_root_length(average_root_len, std_root_len, min_root_len):
    """
    Sample a root sequence length from a normal distribution.
    Mean/std are based on empirical protein length statistics:
    "Protein length distribution is remarkably uniform across the tree of life"
    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02973-2#MOESM2
    """
    root_len = int(np.random.normal(average_root_len, std_root_len))
    return max(root_len, min_root_len)


def sample_indel_params(args):
    root_len = sample_root_length(args.average_root_len, args.std_root_len, args.min_root_len)
    deletion_rate = round(random.uniform(args.min_deletion_rate, args.max_deletion_rate), 4)
    insertion_rate = round(random.uniform(args.min_insertion_rate, args.max_insertion_rate), 4)
    deletion_a_val = round(random.uniform(args.min_deletion_a_val, args.max_deletion_a_val), 4)
    insertion_a_val = round(random.uniform(args.min_insertion_a_val, args.max_insertion_a_val), 4)
    return root_len, deletion_rate, insertion_rate, deletion_a_val, insertion_a_val


def run_simulation(tree, root_len, deletion_rate, insertion_rate, deletion_a_val, insertion_a_val, save_root_sequence):
    """Run one Sailfish simulation and return {node_name: sequence_with_gaps}."""
    sim_protocol = SimProtocol(
        tree,
        root_seq_size=root_len,
        deletion_rate=deletion_rate,
        insertion_rate=insertion_rate,
        deletion_dist=ZipfDistribution(deletion_a_val, 50),
        insertion_dist=ZipfDistribution(insertion_a_val, 50),
        seed=random.randint(0, int(10e7)),
    )

    rep_model = ReplacementModelSpec(
        model=MODEL_CODES.WAG,
        alphabet=ALPHABET_CODES.PROTEIN,
        site_rate_model=SiteRateModelSpec(gamma_alpha=1.0, gamma_categories=4),
    )

    simulator = Simulator(sim_protocol, replacement_model=rep_model)
    if save_root_sequence:
        simulator.save_root_sequence()

    msa = simulator()
    msa_dict = {}
    for i in range(msa.get_num_sequences()):
        name_line, seq = msa.get_msa_row(i).split("\n", 1)
        msa_dict[name_line.replace(">", "").strip()] = seq.strip()
    return msa_dict


def build_source_and_target(data_type, tree, msa_dict, leaf_names, source_delimiter):
    """Return (source_seq, target_seq, special_token, is_tokenizer_source, is_tokenizer_target)."""
    sequences_list = [msa_dict[name] for name in leaf_names]
    source_seq = source_delimiter.join(sequences_list).replace("-", "")

    if data_type == MSA:
        # column-major (transposed) representation of the alignment
        target_seq = "".join(seq[i] for i in range(len(sequences_list[0])) for seq in sequences_list)
        special_token = SPECIAL_MSA
        is_tokenizer_source, is_tokenizer_target = True, True

    elif data_type == TREE:
        target_seq = tree
        special_token = SPECIAL_TREE
        is_tokenizer_source, is_tokenizer_target = True, False

    else:  # ANCESTRAL_SEQUENCE
        target_seq = msa_dict[ANCESTRAL_NODE_NAME].replace("-", "")
        special_token = SPECIAL_ANCESTRAL
        is_tokenizer_source, is_tokenizer_target = True, True

    return source_seq, target_seq, special_token, is_tokenizer_source, is_tokenizer_target


def main():
    args = parse_arguments()
    os.makedirs(args.output_path_folder, exist_ok=True)

    output_file_path = os.path.join(args.output_path_folder, args.output_file_name)
    source_tokenize_path = (
        os.path.join(args.output_path_folder, args.output_source_tokenize_name)
        if args.output_source_tokenize_name else None
    )
    target_tokenize_path = (
        os.path.join(args.output_path_folder, args.output_target_tokenize_name)
        if args.output_target_tokenize_name else None
    )

    leaf_names = [f"l{i}" for i in range(args.num_species)]

    train_writer = BufferedWriter(output_file_path, args.write_every)
    source_writer = BufferedWriter(source_tokenize_path, args.write_every)
    target_writer = BufferedWriter(target_tokenize_path, args.write_every)

    for i in range(args.num_of_simulations):
        if i % 1000 == 0 and i:
            print(f"simulated {i}\t|\t{args.num_of_simulations}")

        tree = build_tree(args, leaf_names)
        root_len, deletion_rate, insertion_rate, deletion_a_val, insertion_a_val = sample_indel_params(args)

        msa_dict = run_simulation(
            tree, root_len, deletion_rate, insertion_rate, deletion_a_val, insertion_a_val,
            save_root_sequence=(args.data_type == ANCESTRAL_SEQUENCE),
        )

        source_seq, target_seq, special_token, is_tokenizer_source, is_tokenizer_target = build_source_and_target(
            args.data_type, tree, msa_dict, leaf_names, args.source_delimiter
        )

        train_writer.add(source_seq + special_token + target_seq)
        if is_tokenizer_source:
            source_writer.add(source_seq)
        if is_tokenizer_target:
            target_writer.add(target_seq)

    train_writer.flush()
    source_writer.flush()
    target_writer.flush()

    special_tokens_path = os.path.join(args.output_path_folder, args.special_tokens_file_name)
    with open(special_tokens_path, "w") as f:
        for token in build_special_tokens_list(leaf_names):
            f.write(f"{token}\n")


if __name__ == "__main__":
    main()
