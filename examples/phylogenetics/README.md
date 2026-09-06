# Phylogenetics Applications — Data Generation Examples

This folder contains two example scripts that use Sailfish to generate
training/evaluation data for common phylogenetics tasks (ancestral sequence
reconstruction, tree inference, and multiple sequence alignment).

---

## `generate_sequences_and_params.py`

Generates pairs of (unaligned sequences, simulation parameters). For each
simulation, a random tree, root length, indel model/rate, and substitution
model are sampled, and Sailfish simulates the resulting sequences. This is
useful for training or testing models that infer simulation parameters
directly from unaligned sequence data.

**Output format** (one line per simulation):

```
<seq_1>|<seq_2>|...|<seq_n>$sub_model: WAG, deletion model: Geo, deletion rate: 0.03, deletion parameter: 0.42, insertion model: Poi, insertion rate: 0.02, insertion parameter: 12.5
```

**Example usage:**

```bash
python generate_sequences_and_params.py \
    --output_path_folder ./data \
    --output_file_name train.txt \
    --num_of_simulations 10000 \
    --num_species 10
```

**Key arguments:**

| Argument | Description |
|---|---|
| `--num_species` | Number of leaves in the simulated tree |
| `--num_of_simulations` | Number of simulations to run |
| `--min/max_branch_length` | Range for random branch lengths |
| `--min/max_deletion_rate`, `--min/max_insertion_rate` | Ranges for indel rates |
| `--average_root_len`, `--std_root_len`, `--min_root_len` | Root sequence length sampling |
| `--source_delimiter` | Delimiter between sequences in the output |
| `--write_every` | How often to flush output to disk |

---

## `generate_msa_tree_ancestral_data.py`

Generates training data for one of three phylogenetics tasks, selected via
`--data_type`:

- **`MSA`** — reconstruct the full multiple sequence alignment (column-major) from unaligned sequences.
- **`ANCESTRAL_SEQUENCE`** — reconstruct the root (ancestral) sequence from unaligned sequences.
- **`TREE`** — infer the guide tree (newick format) from unaligned sequences.

Each output line has the form:

```
<source sequences><special_token><target>
```

where `<special_token>` marks the task (`<ALIGN>`, `<RECONSTRUCT>`, or `<INFER>`
respectively), so a single tokenizer/model can be trained across all three
tasks if desired.

Optionally, the script can also write separate source-only and target-only
files (`--output_source_tokenize_name` / `--output_target_tokenize_name`),
useful for training a tokenizer, plus a `special_tokens.txt` file listing all
leaf names, newick characters, and task tokens.

**Example usage:**

```bash
python generate_msa_tree_ancestral_data.py \
    --data_type ANCESTRAL_SEQUENCE \
    --output_path_folder ./data \
    --output_file_name ancestral_train.txt \
    --output_source_tokenize_name source.txt \
    --output_target_tokenize_name target.txt \
    --num_of_simulations 10000 \
    --num_species 8
```

**Key arguments:**

| Argument | Description |
|---|---|
| `--data_type` | `MSA`, `ANCESTRAL_SEQUENCE`, or `TREE` |
| `--fixed_tree` | Use a fixed 8-species topology instead of a random tree |
| `--num_species` | Number of leaves in the simulated tree |
| `--num_of_simulations` | Number of simulations to run |
| `--min/max_deletion_a_val`, `--min/max_insertion_a_val` | Range for the Zipf distribution shape parameter |
| `--min/max_deletion_rate`, `--min/max_insertion_rate` | Ranges for indel rates |
| `--average_root_len`, `--std_root_len`, `--min_root_len` | Root sequence length sampling |
| `--special_tokens_file_name` | Output file listing special tokens |

---

## `msa_length_distribution.ipynb`

A short notebook demonstrating that, since Sailfish is a Python package, you
can simulate entirely in memory — no need to write MSAs to disk and parse
them back if you only care about a derived statistic.

It fixes a tree, root length, and indel model, runs many replicates, and
plots the distribution of resulting MSA lengths — keeping only the length of
each replicate rather than the full alignment.

Useful as a template for any workflow that needs summary statistics across
many simulations (e.g. exploring how indel parameters affect alignment
length) without the overhead of disk I/O.

---

Both scripts require `msasim`, `ete3`, and `numpy`:

```bash
pip install msasim ete3 numpy
```

The notebook additionally requires `matplotlib`:

```bash
pip install matplotlib
```