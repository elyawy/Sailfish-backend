We provide additional usage examples here. The following examples demonstrate more advanced simulation scenarios using the `msasim` package: partitioned simulations, mixture models, and clade-specific models.

## Partition models

Use `Partitions` to simulate several independent genes/regions (each with its own tree, indel model, and substitution model) and merge them into a single concatenated alignment. Taxa missing from a given partition are gap-filled, so partitions don't need to share the same taxon set.

```python
from msasim import SimProtocol, ReplacementModelSpec, MODEL_CODES, ALPHABET_CODES
from msasim.advanced import Partitions, PartitionModel

# gene 1: fast-evolving, present in 3 taxa
gene1 = PartitionModel(
    name="gene1",
    replacement_model=ReplacementModelSpec(model=MODEL_CODES.GTR, alphabet=ALPHABET_CODES.DNA,
                                            model_parameters=[0.25, 0.25, 0.25, 0.25, 1.0, 4.0, 1.0, 1.0, 4.0, 1.0]),
    indel_model=SimProtocol("(A:0.4,(B:0.3,C:0.3):0.1);", root_seq_size=300, seed=1),
)

# gene 2: slower, present in a different (partially overlapping) taxon set
gene2 = PartitionModel(
    name="gene2",
    replacement_model=ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA),
    indel_model=SimProtocol("(A:0.1,D:0.1);", root_seq_size=450, seed=2),
)

partitions = Partitions()
partitions.add_partition(gene1)
partitions.add_partition(gene2)
partitions.simulate(output_path="concatenated.fasta")
# → taxa A, B, C, D are all written; B/C are gapped across gene2's columns, D is gapped across gene1's columns
```

## Mixture models

Use `MixtureModel` + `Mixture` to simulate a single alignment where each site is drawn from one of several substitution models at random, according to given weights. All models must share the same alphabet. This is useful for approximating rate/model heterogeneity that isn't captured by a single model (e.g. a mix of conserved and fast-evolving site classes).

```python
from msasim import SimProtocol, ReplacementModelSpec, MODEL_CODES, ALPHABET_CODES
from msasim.advanced import MixtureModel, Mixture

conserved = ReplacementModelSpec(model=MODEL_CODES.WAG, alphabet=ALPHABET_CODES.PROTEIN)
fast = ReplacementModelSpec(model=MODEL_CODES.JONES, alphabet=ALPHABET_CODES.PROTEIN)

protocol = SimProtocol("(A:0.2,(B:0.1,C:0.3):0.05);", root_seq_size=500, seed=7)

mixture_model = MixtureModel(
    name="conserved_fast_mix",
    replacement_models=[conserved, fast],
    model_weights=[0.8, 0.2],  # 80% of sites conserved, 20% fast-evolving
    indel_model=protocol,
)

msa = Mixture(mixture_model).simulate()
msa.write_msa("mixture_output.fasta")
```

Note: the alignment's length and gap structure come entirely from the indel model in `MixtureModel.indel_model` — only the substitutions at each site are drawn per-model.

## Clade-specific models

There's currently no dedicated class for clade-specific simulation (unlike `Partitions`/`Mixture` above), but it can be composed from the same building blocks: simulate indels once for the whole tree, simulate substitutions for the whole tree under a base model, then re-simulate substitutions for just one clade under a different model — starting from that clade's own ancestral sequence — and splice the result back into the full alignment. Because the indel history is only simulated once, the result is a single, consistent alignment.

```python
from msasim import SimProtocol, Simulator, ReplacementModelSpec, MODEL_CODES, ALPHABET_CODES

TREE = "(A:0.1,(B:0.1,C:0.1):0.1);"  # clade (B,C) will get a different substitution model
CLADE_TREE = "(B:0.1,C:0.1);"        # same topology/branch lengths as the clade within TREE
CLADE_MRCA_NAME = "N2"               # ete3/Sailfish's default label for the (B,C) ancestor
ROOT_LEN = 300

# 1. Simulate indels + substitutions (model A) for the whole tree in one pass.
#    save_all_nodes_sequences() also gives us the clade's ancestral sequence.
protocol = SimProtocol(TREE, root_seq_size=ROOT_LEN, deletion_rate=0.05, insertion_rate=0.05, seed=1)
model_a = ReplacementModelSpec(model=MODEL_CODES.WAG, alphabet=ALPHABET_CODES.PROTEIN)
sim_a = Simulator(protocol, replacement_model=model_a)
sim_a.save_all_nodes_sequences()
msa_a = sim_a()

rows = {}
for i in range(msa_a.get_num_sequences()):
    name, seq = msa_a.get_msa_row(i).split("\n", 1)
    rows[name.lstrip(">")] = seq  # includes leaves and internal nodes, e.g. "N2"

clade_mrca_gapped = rows[CLADE_MRCA_NAME]
clade_mrca_seq = clade_mrca_gapped.replace("-", "")  # strip alignment gaps to get a plain sequence

# 2. Re-simulate substitutions only for the clade, starting from its ancestral sequence,
#    under a different model. No further indels — the clade keeps the gap pattern from step 1.
clade_protocol = SimProtocol(CLADE_TREE, root_seq_size=len(clade_mrca_seq),
                              deletion_rate=0.0, insertion_rate=0.0, seed=2)
model_b = ReplacementModelSpec(model=MODEL_CODES.JONES, alphabet=ALPHABET_CODES.PROTEIN)
sim_b = Simulator(clade_protocol, replacement_model=model_b)
sim_b.set_root_sequence(clade_mrca_seq)
msa_b = sim_b()

clade_rows = {}
for i in range(msa_b.get_num_sequences()):
    name, seq = msa_b.get_msa_row(i).split("\n", 1)
    clade_rows[name.lstrip(">")] = seq


def reinsert_gaps(template_with_gaps, new_ungapped_chars):
    """Re-insert '-' at the same positions as in template_with_gaps, filling the rest from new_ungapped_chars."""
    chars = iter(new_ungapped_chars)
    return "".join("-" if c == "-" else next(chars) for c in template_with_gaps)


# 3. Build the final single alignment: model-A sequences outside the clade,
#    model-B sequences (re-gapped to match) inside it.
final_msa = {"A": rows["A"]}
for leaf in ("B", "C"):
    final_msa[leaf] = reinsert_gaps(rows[leaf], clade_rows[leaf])

with open("clade_specific.fasta", "w") as f:
    for name, seq in final_msa.items():
        f.write(f">{name}\n{seq}\n")
```

This gives one alignment, with a single shared indel history, where clade `(B, C)` evolved under `model_b` while the rest of the tree evolved under `model_a`. The recipe generalizes to nested or multiple clades by repeating step 2 for each clade's MRCA name and model, and to deeper subtrees by extracting the correct clade topology/branch lengths from the original tree (e.g. with `ete3`) instead of hardcoding `CLADE_TREE`.