> NOTICE: SAILFISH IS IN EARLY DEVELOPMENT AND SUBJECT TO CHANGE. API AND FEATURES MAY BE UNSTABLE. IT IS PROVIDED AS-IS WITHOUT ANY WARRANTY. USE AT YOUR OWN RISK. PLEASE REPORT ISSUES AND SUGGESTIONS ON GITHUB.

# Sailfish

Sailfish is a high-performance multiple sequence alignment (MSA) simulator written in C++ with a Python API. It enables rapid generation of large-scale simulated datasets with support for indels, substitutions, and realistic evolutionary models.

## Features

- High-performance C++ engine with ergonomic Python interface
- Support for DNA, protein, and codon sequence evolution
- Flexible indel modeling with multiple length distributions (Zipf, Geometric, Poisson, Custom)
- 26+ substitution models including JTT, WAG, LG, HKY, GTR, and more
- Gamma rate heterogeneity, invariant sites, and site rate correlation
- Per-branch parameter specification for heterogeneous models
- Low-memory mode for large-scale simulations (1M+ sequences)
- Reproducible simulations with explicit seed control

## Installation

```bash
pip install msasim
```

Requirements: Python >= 3.7

## Web app

[Try Sailfish without installing anything](https://elyawy.github.io/SailfishWeb/) — runs entirely in your browser via WASM.


## Quick Start

### Full simulation (indels + substitutions)

```python
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

### Indels only (no substitutions)

Omit the `replacement_model` argument:

```python
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

### Substitutions only (no indels)

Set both rates to zero:

```python
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

### Gamma rate heterogeneity

```python
from msasim import SimProtocol, Simulator, ReplacementModelSpec, SiteRateModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES

protocol = SimProtocol(tree="(A:0.5,B:0.5);", root_seq_size=100, seed=42)

rate_model = SiteRateModelSpec(gamma_alpha=1.0, gamma_categories=4)
rep_model = ReplacementModelSpec(
    model=MODEL_CODES.WAG,
    alphabet=ALPHABET_CODES.PROTEIN,
    site_rate_model=rate_model,
)
simulator = Simulator(protocol, replacement_model=rep_model)
msa = simulator()
```

### Batch simulations

```python
from msasim import SimProtocol, Simulator, ReplacementModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES

protocol = SimProtocol(tree="tree.nwk", root_seq_size=500, seed=42)
rep_model = ReplacementModelSpec(model=MODEL_CODES.LG, alphabet=ALPHABET_CODES.PROTEIN)
simulator = Simulator(protocol, replacement_model=rep_model)

# Internal RNG advances deterministically across calls
for i in range(100):
    msa = simulator()
    msa.write_msa(f"replicate_{i:04d}.fasta")
```

### Low-memory mode for large simulations

```python
import pathlib
from msasim import SimProtocol, Simulator, ReplacementModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES

protocol = SimProtocol(tree="large_tree.nwk", root_seq_size=10000, seed=42)
rep_model = ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA)
simulator = Simulator(protocol, replacement_model=rep_model)

simulator.simulate(output_path=pathlib.Path("large_alignment.fasta"))
# → writes large_alignment_replicate_1.fasta
```

## Documentation

For the complete API reference see [API_REFERENCE.md](API_REFERENCE.md).

## Core Concepts

### Simulation type is inferred from the model

There is no `simulation_type` argument. The alphabet (`DNA`, `PROTEIN`, `CODON`) is taken from the `ReplacementModelSpec`. Passing no `replacement_model` runs an indels-only simulation.

### Indel length distributions

- `ZipfDistribution(a, truncation)` — power-law, typical for biological indels
- `GeometricDistribution(p, truncation)` — exponentially decreasing
- `PoissonDistribution(lambda, truncation)` — Poisson-based
- `CustomDistribution(probs)` — user-defined probability vector

### Available substitution models

**DNA:** `NUCJC`, `HKY`, `GTR`, `TAMURA92`  
**Protein:** `WAG`, `LG`, `JONES` (JTT), `DAYHOFF`, `MTREV24`, `CPREV45`, `HIVB`, `HIVW`, and more  
**Codon:** `EMPIRICODON`, `WYANG`

See [API_REFERENCE.md](API_REFERENCE.md) for the full list.


## Project Goals

- Ease of use: Simple, intuitive Python API
- Speed: High-performance C++ implementation
- Modularity: Flexible configuration of all evolutionary parameters

## Contributing

Bug reports and feature requests are welcome via GitHub issues.

## Citation

If you use Sailfish in your research, please cite:

[Citation information to be added]

## License

Academic Free License v3.0. See [LICENSE](LICENSE) for details.
