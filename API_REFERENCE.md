# msasim API Reference

Complete reference for the msasim Python API.

## Table of Contents

- [Installation](#installation)
- [Imports](#imports)
- [SimProtocol](#simprotocol)
- [Simulator](#simulator)
- [ReplacementModelSpec](#replacementmodelspec)
- [SiteRateModelSpec](#siteratemodespec)
- [Distributions](#distributions)
- [Msa](#msa)
- [Substitution Models](#substitution-models)
- [Advanced: Partitions, Mixtures, and Custom Models](#advanced-partitions-mixtures-and-custom-models)

---

## Installation

```bash
pip install msasim
```

---

## Imports

```python
from msasim import (
    SimProtocol,
    Simulator,
    Msa,
    ReplacementModelSpec,
    SiteRateModelSpec,
    MODEL_CODES,
    ALPHABET_CODES,
    ZipfDistribution,
    GeometricDistribution,
    PoissonDistribution,
    CustomDistribution,
)
```

> **Note:** `from msasim import sailfish` is deprecated. Import directly from `msasim`.

---

## SimProtocol

Configuration object for a simulation: tree, indel rates, distributions, seed.

### Constructor

```python
SimProtocol(
    tree: Union[Tree, str],          # newick string, file path, or Tree object
    root_seq_size: int = 100,
    deletion_rate: float = 0.05,
    insertion_rate: float = 0.05,
    deletion_dist: Distribution = CustomDistribution([0.5, 0.3, 0.2]),
    insertion_dist: Distribution = CustomDistribution([0.5, 0.3, 0.2]),
    minimum_seq_size: int = 100,
    seed: Optional[int] = None,      # defaults to time_ns() if omitted
)
```

### Key setter methods

```python
protocol.set_sequence_size(n: int) -> None
protocol.set_seed(seed: int) -> None
protocol.set_insertion_rates(rate: float) -> None          # uniform across branches
protocol.set_insertion_rates(rates: List[float]) -> None   # per-branch
protocol.set_deletion_rates(rate: float) -> None
protocol.set_deletion_rates(rates: List[float]) -> None
protocol.set_insertion_length_distributions(dist: Distribution) -> None
protocol.set_deletion_length_distributions(dist: Distribution) -> None
protocol.set_min_sequence_size(n: int) -> None
protocol.set_cache_branch_probs(cache: bool) -> None  # precompute branch probs; speeds up multi-replicate runs
```

### Example

```python
from msasim import SimProtocol
from msasim.distributions import ZipfDistribution

protocol = SimProtocol(
    tree="((A:0.1,B:0.1):0.2,C:0.3);",
    root_seq_size=500,
    deletion_rate=0.05,
    insertion_rate=0.02,
    deletion_dist=ZipfDistribution(1.7, 50),
    insertion_dist=ZipfDistribution(1.7, 50),
    seed=42,
)
```

---

## Simulator

Main simulation engine.

### Constructor

```python
Simulator(
    simProtocol: SimProtocol,
    replacement_model: Optional[ReplacementModelSpec] = None,
)
```

- Pass a `ReplacementModelSpec` to run substitutions. The alphabet (`DNA`/`PROTEIN`/`CODON`) is taken from it.
- Omit `replacement_model` (or pass `None`) to run an indels-only simulation.
- Whether indels are generated is determined by whether the rates in the protocol are non-zero — there is no separate flag.

### Running a simulation

```python
# Single run — returns an Msa
msa = simulator()

# Multiple replicates — returns List[Msa]
msas = simulator.simulate(times=100)

# Low-memory mode — writes to disk, returns None
# Files are named {stem}_replicate_1.fasta, {stem}_replicate_2.fasta, etc.
simulator.simulate(output_path=pathlib.Path("output.fasta"))
# → writes output_replicate_1.fasta
```

### Substitution model

To swap the model after construction (must keep the same alphabet):

```python
simulator.set_replacement_model(replacement_model: ReplacementModelSpec) -> None
```

### Sequence saving

By default only leaf sequences are saved.

```python
simulator.save_root_sequence() -> None
simulator.save_all_nodes_sequences() -> None
simulator.save_leaves_sequences() -> None
simulator.get_sequences_to_save() -> List[bool]
```

### Per-site rates

```python
simulator.save_rates(True)           # must call before simulate
msa = simulator()
rates = simulator.get_rates()        # List[float], one per root site
categories = simulator.get_rate_categories()  # List[int], rate category index per site
```

### Root sequence

```python
simulator.save_root_sequence()
simulator.set_root_sequence("ACGT...")   # length must match root_seq_size
```

### Examples

```python
# Indels only
from msasim import SimProtocol, Simulator
from msasim.distributions import ZipfDistribution

protocol = SimProtocol(
    tree="(A:0.5,B:0.5);",
    root_seq_size=100,
    deletion_rate=0.05,
    insertion_rate=0.02,
    deletion_dist=ZipfDistribution(1.7, 50),
    insertion_dist=ZipfDistribution(1.7, 50),
    seed=42,
)
simulator = Simulator(protocol)
msa = simulator()
```

```python
# Substitutions only (no indels)
from msasim import SimProtocol, Simulator, ReplacementModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES

protocol = SimProtocol(
    tree="(A:0.5,B:0.5);",
    root_seq_size=100,
    insertion_rate=0.0,
    deletion_rate=0.0,
    seed=42,
)
rep_model = ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA)
simulator = Simulator(protocol, replacement_model=rep_model)
msa = simulator()
```

```python
# Full simulation (indels + substitutions)
from msasim import SimProtocol, Simulator, ReplacementModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES
from msasim.distributions import ZipfDistribution

protocol = SimProtocol(
    tree="(A:0.5,B:0.5);",
    root_seq_size=100,
    deletion_rate=0.05,
    insertion_rate=0.02,
    deletion_dist=ZipfDistribution(1.7, 50),
    insertion_dist=ZipfDistribution(1.7, 50),
    seed=42,
)
rep_model = ReplacementModelSpec(model=MODEL_CODES.WAG, alphabet=ALPHABET_CODES.PROTEIN)
simulator = Simulator(protocol, replacement_model=rep_model)
msa = simulator()
msa.write_msa("output.fasta")
```

---

## ReplacementModelSpec

Specifies the substitution model and (optionally) the site rate model.

```python
ReplacementModelSpec(
    model: MODEL_CODES,
    alphabet: ALPHABET_CODES,
    model_parameters: Optional[List[float]] = None,  # required for HKY, GTR, etc.
    amino_model_file: Optional[pathlib.Path] = None, # for REVERSIBLE custom protein model
    site_rate_model: SiteRateModelSpec = SiteRateModelSpec(),
)
```

### Examples

```python
from msasim import ReplacementModelSpec, SiteRateModelSpec, MODEL_CODES, ALPHABET_CODES

# JC — no parameters
rep_model = ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA)

# HKY — requires frequencies + kappa
rep_model = ReplacementModelSpec(
    model=MODEL_CODES.HKY,
    alphabet=ALPHABET_CODES.DNA,
    model_parameters=[0.25, 0.25, 0.25, 0.25, 2.0],  # freqs A/C/G/T, kappa
)

# WAG protein with gamma rate heterogeneity
rate_model = SiteRateModelSpec(gamma_alpha=1.0, gamma_categories=4)
rep_model = ReplacementModelSpec(
    model=MODEL_CODES.WAG,
    alphabet=ALPHABET_CODES.PROTEIN,
    site_rate_model=rate_model,
)
```

---

## SiteRateModelSpec

Controls per-site rate variation. Used inside `ReplacementModelSpec`.

```python
SiteRateModelSpec(
    gamma_alpha: float = 1.0,
    gamma_categories: int = 1,
    invariant_proportion: float = 0.0,
    site_rate_correlation: float = 0.0,
    free_rates: Optional[List[float]] = None,
    free_rate_weights: Optional[List[float]] = None,
)
```

- `gamma_alpha`: shape parameter of the discretized gamma distribution; lower values = more rate variation.
- `gamma_categories`: number of discrete rate categories.
- `invariant_proportion`: fraction of sites that are completely invariant (evolve at rate zero). Must be in `[0, 1)`.
- `site_rate_correlation`: autocorrelation between adjacent sites' rate categories. Must be in `[0, 1)`; requires `gamma_categories > 1`.
- `free_rates` / `free_rate_weights`: explicit rate values and their weights (must sum to 1.0). Mutually exclusive with the gamma model; if provided, free rates take precedence.

### Example

```python
from msasim import SiteRateModelSpec

rate_model = SiteRateModelSpec(
    gamma_alpha=0.5,
    gamma_categories=8,
    site_rate_correlation=0.99,
)
```

---

## Distributions

All distributions take a `truncation` argument limiting the maximum indel length.

```python
ZipfDistribution(a: float, truncation: int = 150)
GeometricDistribution(p: float, truncation: int = 150)
PoissonDistribution(lam: float, truncation: int = 150)
CustomDistribution(probs: List[float])   # probs[i] = P(length = i+1), must sum to 1
```

---

## Msa

The result object returned by `simulator()`.

```python
msa.get_length() -> int                  # alignment length (including gaps)
msa.get_num_sequences() -> int
msa.get_msa_row(i: int) -> str           # ">name\nsequence"
msa.get_sparse_msa() -> List            # compact block representation
msa.get_per_site_rate_categories() -> List[int]
msa.write_msa(path: str) -> None        # write FASTA to file
msa.print_msa() -> None
```

---

## Substitution Models

### DNA (`ALPHABET_CODES.DNA`)

| `MODEL_CODES` | Parameters (`model_parameters`) |
|---|---|
| `NUCJC` | none |
| `HKY` | `[freq_A, freq_C, freq_G, freq_T, kappa]` |
| `GTR` | `[freq_A, freq_C, freq_G, freq_T, rate_AC, rate_AG, rate_AT, rate_CG, rate_CT, rate_GT]` |
| `TAMURA92` | — |

### Protein (`ALPHABET_CODES.PROTEIN`)

`WAG`, `LG`, `JONES` (JTT), `DAYHOFF`, `MTREV24`, `CPREV45`, `HIVB`, `HIVW`, `AAJC`,
`EX_BURIED`, `EX_EXPOSED`, `EHO_EXTENDED`, `EHO_HELIX`, `EHO_OTHER`,
`EX_EHO_BUR_EXT`, `EX_EHO_BUR_HEL`, `EX_EHO_BUR_OTH`, `EX_EHO_EXP_EXT`, `EX_EHO_EXP_HEL`, `EX_EHO_EXP_OTH`

### Codon (`ALPHABET_CODES.CODON`)

| `MODEL_CODES` | Parameters |
|---|---|
| `EMPIRICODON` | none |
| `WYANG` | `[selection_factor, kappa]` |

---

## Advanced: Partitions, Mixtures, and Custom Models

### Custom and Non-reversible Models

`REVERSIBLE` and `NONREVERSIBLE` are available for all alphabet types (`DNA`, `PROTEIN`, `CODON`) but require special configuration.

**`REVERSIBLE`** — load a custom time-reversible model from a file:

```python
rep_model = ReplacementModelSpec(
    model=MODEL_CODES.REVERSIBLE,
    alphabet=ALPHABET_CODES.PROTEIN,
    amino_model_file=pathlib.Path("my_model.txt"),
)
```

**`NONREVERSIBLE`** — provide the full instantaneous rate matrix Q and equilibrium frequencies directly in `model_parameters`. The expected length is `alphabet_size * alphabet_size + alphabet_size` (Q matrix entries row-by-row, then frequencies):

```python
# Example for DNA (alphabet_size=4): 4*4 + 4 = 20 values
rep_model = ReplacementModelSpec(
    model=MODEL_CODES.NONREVERSIBLE,
    alphabet=ALPHABET_CODES.DNA,
    model_parameters=[
        # Q matrix (4x4, row by row)
        0.0,  0.1,  0.2,  0.3,
        0.4,  0.0,  0.1,  0.2,
        0.1,  0.3,  0.0,  0.1,
        0.2,  0.1,  0.3,  0.0,
        # equilibrium frequencies
        0.25, 0.25, 0.25, 0.25,
    ],
)
```

### Partitions and Mixtures

#### Partitions

Simulate multiple partitions over different taxa sets and merge into a single FASTA.

```python
from msasim import SimProtocol, ReplacementModelSpec, MODEL_CODES, ALPHABET_CODES
from msasim.advanced import Partitions, PartitionModel

partition1 = PartitionModel(
    name="p1",
    replacement_model=ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA),
    indel_model=SimProtocol("(A:0.1,B:0.1);", root_seq_size=200, seed=1),
)
partition2 = PartitionModel(
    name="p2",
    replacement_model=ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA),
    indel_model=SimProtocol("(A:0.1,C:0.1);", root_seq_size=150, seed=2),
)

partitions = Partitions()
partitions.add_partition(partition1)
partitions.add_partition(partition2)
partitions.simulate(output_path="merged.fasta")
```

#### Mixtures

Simulate a mixture of substitution models where each site is drawn from one model at random.

```python
from msasim import SimProtocol, ReplacementModelSpec, MODEL_CODES, ALPHABET_CODES
from msasim.advanced import MixtureModel, Mixture

models = [ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA) for _ in range(3)]
weights = [0.3, 0.3, 0.4]

protocol = SimProtocol("(A:0.1,B:0.1,C:0.1);", root_seq_size=100, seed=1)
mixture_model = MixtureModel(
    name="my_mixture",
    replacement_models=models,
    model_weights=weights,
    indel_model=protocol,
)
msa = Mixture(mixture_model).simulate()
```