
# Sailfish

Sailfish is a high-performance multiple sequence alignment (MSA) simulator written in C++ with a Python API. It enables rapid generation of large-scale simulated datasets with support for indels, substitutions, and realistic evolutionary models.

## Features

- High-performance C++ engine with ergonomic Python interface
- Support for DNA, protein, and codon sequence evolution
- Flexible indel modeling with multiple length distributions (Zipf, Geometric, Poisson, Custom)
- 26+ substitution models including JTT, WAG, LG, HKY, GTR, and more
- Gamma rate heterogeneity, invariant sites, and site rate correlation
- Per-branch parameter specification for heterogeneous models
- Low-memory mode for large-scale simulations (1M+ sites)
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

# Define a simulation protocol with a Newick tree, root sequence size, and indel parameters
protocol = SimProtocol(
    tree="(A:0.5,B:0.5);", # Newick tree string or path to a Newick file
    root_seq_size=100, # number of sites in the root sequence
    deletion_rate=0.05, # rate of deletions per substitutions
    insertion_rate=0.02, # rate of insertions per substitutions
    deletion_dist=ZipfDistribution(1.7, 50), # indel length distribution for deletions (Zipf distribution)
    insertion_dist=ZipfDistribution(1.7, 50), # indel length distribution for insertions (Zipf distribution)
)

# Define a substitution model (WAG for proteins in this case)
rep_model = ReplacementModelSpec(model=MODEL_CODES.WAG, alphabet=ALPHABET_CODES.PROTEIN)

# Create a simulator instance with the protocol and substitution model
simulator = Simulator(protocol, replacement_model=rep_model)

# Run the simulation to generate a multiple sequence alignment (MSA)
msa = simulator()
msa.write_msa("output.fasta")
```

## Documentation

For the complete API reference see [API_REFERENCE.md](API_REFERENCE.md).
For usage examples see [examples/README.md](examples/README.md).

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
