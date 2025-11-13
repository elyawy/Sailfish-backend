# Sailfish

Sailfish is a high-performance multiple sequence alignment (MSA) simulator written in C++ with a Python API. It enables rapid generation of large-scale simulated datasets with support for indels, substitutions, and realistic evolutionary models.

## Features

- High-performance C++ engine with ergonomic Python interface
- Support for both DNA and protein sequence evolution
- Flexible indel modeling with multiple length distributions (Zipf, Geometric, Poisson, Custom)
- 26+ substitution models including JTT, WAG, LG, HKY, GTR, and more
- Gamma rate heterogeneity and invariant sites
- Per-branch parameter specification for heterogeneous models
- Low-memory mode for large-scale simulations (1M+ sequences)
- Reproducible simulations with explicit seed control

## Installation

```bash
pip install msasim
```

Requirements: Python >= 3.6

## Quick Start

### Basic Example with Indels and Substitutions

```python
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES, ZipfDistribution

# Configure simulation protocol
sim_protocol = sim.SimProtocol(
    tree="(A:0.5,B:0.5);",
    root_seq_size=100,
    deletion_rate=0.01,
    insertion_rate=0.01,
    deletion_dist=ZipfDistribution(1.7, 50),
    insertion_dist=ZipfDistribution(1.7, 50),
    seed=42
)

# Create simulator
simulation = sim.Simulator(sim_protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)

# Configure substitution model with gamma rate heterogeneity
simulation.set_replacement_model(
    model=MODEL_CODES.WAG,
    gamma_parameters_alpha=1.0,
    gamma_parameters_categories=4
)

# Run simulation
msa = simulation()

# Output results
msa.write_msa("output.fasta")
msa.print_msa()
```

### Substitutions-Only Simulation

```python
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES

# No indels configured
protocol = sim.SimProtocol(
    tree="path/to/tree.nwk",
    root_seq_size=500,
    seed=42
)

simulator = sim.Simulator(protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)
simulator.set_replacement_model(model=MODEL_CODES.LG)

msa = simulator()
msa.write_msa("alignment.fasta")
```

### Batch Simulations

```python
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES

# Initialize once with seed
protocol = sim.SimProtocol(tree="tree.nwk", root_seq_size=500, seed=42)
simulator = sim.Simulator(protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)
simulator.set_replacement_model(model=MODEL_CODES.JTT)

# Generate multiple replicates
# Internal RNG advances automatically for reproducibility
for i in range(100):
    msa = simulator()
    msa.write_msa(f"replicate_{i:04d}.fasta")
```

### Low-Memory Mode for Large Simulations

```python
import pathlib
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES

protocol = sim.SimProtocol(tree="large_tree.nwk", root_seq_size=10000, seed=42)
simulator = sim.Simulator(protocol, simulation_type=sim.SIMULATION_TYPE.DNA)
simulator.set_replacement_model(model=MODEL_CODES.NUCJC)

# Write directly to disk without holding MSA in memory
simulator.simulate_low_memory(pathlib.Path("large_alignment.fasta"))
```

## Documentation

For complete API documentation, including all available models, distributions, and advanced features, see [API_REFERENCE.md](API_REFERENCE.md).

## Core Concepts

### Simulation Types

- `SIMULATION_TYPE.NOSUBS`: Indels only, no substitutions
- `SIMULATION_TYPE.DNA`: DNA sequences with nucleotide models
- `SIMULATION_TYPE.PROTEIN`: Protein sequences with amino acid models

### Available Models

**Nucleotide**: JC, HKY, GTR, Tamura92
**Protein**: WAG, LG, JTT (JONES), Dayhoff, MTREV24, CPREV45, HIV models, and more

See [API_REFERENCE.md#substitution-models](API_REFERENCE.md#substitution-models) for the complete list.

### Indel Length Distributions

- **ZipfDistribution**: Power-law distribution (typical for biological data)
- **GeometricDistribution**: Exponentially decreasing
- **PoissonDistribution**: Poisson-based
- **CustomDistribution**: User-defined probability vector

## Common Use Cases

### 1. Phylogenetic Method Validation

Generate known-truth alignments with specified evolutionary parameters to test inference methods.

### 2. Benchmarking Alignment Tools

Create challenging datasets with varying indel rates and substitution patterns.

### 3. Statistical Power Analysis

Simulate datasets under different evolutionary scenarios to assess method sensitivity.

### 4. Model Comparison

Generate alignments under different substitution models for model selection studies.

## Performance Notes

Typical simulation times on modern hardware:
- 10K sequences × 30K sites: ~10 seconds
- 100K sequences × 30K sites: ~2 minutes
- 1M sequences × 30K sites: ~20 minutes

For simulations with >100K sequences or memory constraints, use `simulate_low_memory()`.

Memory usage estimate: `(num_sequences × alignment_length) / 300,000` MB

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

[License information to be added]