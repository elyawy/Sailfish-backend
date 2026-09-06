"""
Generate (unaligned sequences, simulation parameters) pairs for training/testing.

For each simulation:
  1. A random tree topology is generated and branch lengths are drawn uniformly.
  2. Root sequence length is sampled from a normal distribution based on
     empirical protein length statistics (see reference below).
  3. Indel model, rates and parameters are sampled for insertion and deletion
     independently.
  4. A substitution model is sampled and the simulation is run with Sailfish.
  5. Output format per line:
       <unaligned sequences joined by delimiter><SEP><parameters string>
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
    GeometricDistribution,
    PoissonDistribution,
)

# separates the input sequences from the target parameter string in each output line
SEQUENCE_PARAMS_SEPARATOR = "$"

SUBSTITUTION_MODELS = [
    ("DAYHOFF", MODEL_CODES.DAYHOFF),
    ("LG", MODEL_CODES.LG),
    ("JONES", MODEL_CODES.JONES),
    ("MTREV24", MODEL_CODES.MTREV24),
    ("WAG", MODEL_CODES.WAG),
    ("CPREV45", MODEL_CODES.CPREV45),
    ("HIVB", MODEL_CODES.HIVB),
    ("HIVW", MODEL_CODES.HIVW),
]

# (name, distribution class, (param_min, param_max), truncation_length)
INDEL_LENGTH_DISTRIBUTIONS = [
    ("Zipf", ZipfDistribution, (1.03, 2.78), 150),
    ("Geo", GeometricDistribution, (0.04, 0.66), 150),
    ("Poi", PoissonDistribution, (1.5, 25), 150),
]


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate simulated sequence/parameter training data")

    # files
    parser.add_argument("--output_path_folder", type=str, required=True, help="output folder")
    parser.add_argument("--output_file_name", type=str, required=True, help="output file name")

    # tree
    parser.add_argument("--num_species", type=int, default=10, help="number of species (leaves)")
    parser.add_argument("--min_branch_length", type=float, default=0.00, help="minimum branch length")
    parser.add_argument("--max_branch_length", type=float, default=1.00, help="maximum branch length")

    # simulator settings
    parser.add_argument("--num_of_simulations", type=int, required=True, help="number of simulations to run")
    parser.add_argument("--min_deletion_rate", type=float, default=0.00, help="minimum deletion rate")
    parser.add_argument("--max_deletion_rate", type=float, default=0.05, help="maximum deletion rate")
    parser.add_argument("--min_insertion_rate", type=float, default=0.00, help="minimum insertion rate")
    parser.add_argument("--max_insertion_rate", type=float, default=0.05, help="maximum insertion rate")
    parser.add_argument("--average_root_len", type=int, default=345, help="mean of root sequence length")
    parser.add_argument("--std_root_len", type=int, default=76, help="std of root sequence length")
    parser.add_argument("--min_root_len", type=int, default=25, help="lower bound on root sequence length")

    # output formatting
    parser.add_argument("--source_delimiter", type=str, default="|", help="delimiter between sequences")
    parser.add_argument("--write_every", type=int, default=1000, help="flush output every N simulations")

    return parser.parse_args()


def build_random_tree(leaf_names, min_branch_length, max_branch_length):
    """Build a random topology over the given leaf names with uniform random branch lengths."""
    t = Tree()
    shuffled_names = leaf_names.copy()
    random.shuffle(shuffled_names)
    t.populate(len(leaf_names), names_library=shuffled_names)

    newick = t.write(format=5).replace(":1", "##")
    while "##" in newick:
        branch_length = round(random.uniform(min_branch_length, max_branch_length), 3)
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


def sample_indel_setting(rate_min, rate_max):
    """Pick an indel length distribution, its rate and its parameter value."""
    name, dist_class, param_range, truncation = random.choice(INDEL_LENGTH_DISTRIBUTIONS)
    rate = round(random.uniform(rate_min, rate_max), 4)
    param_value = round(random.uniform(*param_range), 4)
    return name, dist_class, truncation, rate, param_value


def run_single_simulation(args, leaf_names):
    """Run one simulation and return (unaligned sequences dict, parameters description string)."""
    tree = build_random_tree(leaf_names, args.min_branch_length, args.max_branch_length)
    root_len = sample_root_length(args.average_root_len, args.std_root_len, args.min_root_len)

    d_name, d_dist_class, d_truncation, deletion_rate, d_param = sample_indel_setting(
        args.min_deletion_rate, args.max_deletion_rate
    )
    i_name, i_dist_class, i_truncation, insertion_rate, i_param = sample_indel_setting(
        args.min_insertion_rate, args.max_insertion_rate
    )

    sim_protocol = SimProtocol(
        tree,
        root_seq_size=root_len,
        deletion_rate=deletion_rate,
        insertion_rate=insertion_rate,
        deletion_dist=d_dist_class(d_param, d_truncation),
        insertion_dist=i_dist_class(i_param, i_truncation),
        seed=random.randint(0, int(10e7)),
    )

    sub_model_name, sub_model = random.choice(SUBSTITUTION_MODELS)
    rep_model = ReplacementModelSpec(
        model=sub_model,
        alphabet=ALPHABET_CODES.PROTEIN,
        site_rate_model=SiteRateModelSpec(gamma_alpha=1.0, gamma_categories=4),
    )

    simulator = Simulator(sim_protocol, replacement_model=rep_model)
    msa = simulator()

    msa_dict = {}
    for i in range(msa.get_num_sequences()):
        name_line, seq = msa.get_msa_row(i).split("\n", 1)
        msa_dict[name_line.replace(">", "").strip()] = seq.strip()

    params_str = (
        f"sub_model: {sub_model_name}, "
        f"deletion model: {d_name}, deletion rate: {deletion_rate}, deletion parameter: {d_param}, "
        f"insertion model: {i_name}, insertion rate: {insertion_rate}, insertion parameter: {i_param}"
    )
    return msa_dict, params_str


def flush_to_file(output_file_path, lines):
    with open(output_file_path, "a") as f:
        for line in lines:
            f.write(line + "\n")


def main():
    args = parse_arguments()
    os.makedirs(args.output_path_folder, exist_ok=True)
    output_file_path = os.path.join(args.output_path_folder, args.output_file_name)

    leaf_names = [f"l{i}" for i in range(args.num_species)]

    buffered_lines = []
    for i in range(args.num_of_simulations):
        if i % 1000 == 0 and i:
            print(f"simulated {i}\t|\t{args.num_of_simulations}")

        msa_dict, params_str = run_single_simulation(args, leaf_names)

        # unaligned: drop gaps, keep leaf order fixed
        sequences = [msa_dict[name] for name in leaf_names]
        source_seq = args.source_delimiter.join(sequences).replace("-", "")

        buffered_lines.append(source_seq + SEQUENCE_PARAMS_SEPARATOR + params_str)

        if len(buffered_lines) >= args.write_every:
            flush_to_file(output_file_path, buffered_lines)
            buffered_lines = []

    flush_to_file(output_file_path, buffered_lines)


if __name__ == "__main__":
    main()