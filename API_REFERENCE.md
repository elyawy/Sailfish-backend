# Sailfish (msasim) API Reference

Complete reference documentation for the Sailfish Python API. Sailfish is a high-performance multiple sequence alignment simulator with support for indels and substitutions.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Core Concepts](#core-concepts)
- [API Reference](#api-reference)
  - [SimProtocol](#simprotocol)
  - [Simulator](#simulator)
  - [Distributions](#distributions)
  - [Tree](#tree)
  - [Msa](#msa)
- [Substitution Models](#substitution-models)
- [Usage Patterns](#usage-patterns)
- [Advanced Features](#advanced-features)
- [Performance Considerations](#performance-considerations)

## Installation

```bash
pip install msasim
```

Requirements:
- Python >= 3.6
- C++ compiler with C++17 support (for building from source)

## Quick Start

```python
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES, ZipfDistribution

# Create simulation protocol
protocol = sim.SimProtocol(
    tree="(A:0.5,B:0.5);",
    root_seq_size=100,
    deletion_rate=0.01,
    insertion_rate=0.01,
    deletion_dist=ZipfDistribution(1.7, 50),
    insertion_dist=ZipfDistribution(1.7, 50),
    seed=42
)

# Create simulator
simulator = sim.Simulator(protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)

# Configure substitution model
simulator.set_replacement_model(
    model=MODEL_CODES.WAG,
    gamma_parameters_alpha=1.0,
    gamma_parameters_categories=4
)

# Run simulation
msa = simulator()

# Output results
msa.write_msa("output.fasta")
```

## Core Concepts

### Workflow

Sailfish simulations follow a four-step workflow:

1. **Protocol Definition**: Configure tree, indel rates, and distributions
2. **Simulator Creation**: Specify simulation type (DNA/Protein/NoSubs)
3. **Model Configuration**: Set substitution model and parameters
4. **Execution**: Generate MSA(s) and output results

### Simulation Types

Three simulation types are available via `sim.SIMULATION_TYPE`:

- `NOSUBS`: Indels only, no substitutions
- `DNA`: DNA sequences with nucleotide substitution models
- `PROTEIN`: Protein sequences with amino acid substitution models

### Parameter Hierarchy

Parameters can be specified at three levels:

1. **Constructor**: Initial values in `SimProtocol` constructor
2. **Setter methods**: Per-parameter configuration via `set_*` methods
3. **Per-branch**: Different values for each branch (optional)

## API Reference

### SimProtocol

Configuration object that defines simulation parameters.

#### Constructor

```python
SimProtocol(
    tree: Union[Tree, str],
    root_seq_size: int = 100,
    deletion_rate: float = 0.0,
    insertion_rate: float = 0.0,
    deletion_dist: Distribution = ZipfDistribution(1.7, 50),
    insertion_dist: Distribution = ZipfDistribution(1.7, 50),
    minimum_seq_size: int = 100,
    seed: int = 0
)
```

**Parameters:**
- `tree`: Phylogenetic tree in Newick format (string) or file path
- `root_seq_size`: Length of root sequence (amino acids or nucleotides)
- `deletion_rate`: Per-site deletion rate per unit branch length
- `insertion_rate`: Per-site insertion rate per unit branch length
- `deletion_dist`: Length distribution for deletions
- `insertion_dist`: Length distribution for insertions
- `minimum_seq_size`: Minimum sequence length constraint
- `seed`: Random seed for reproducibility

**Example:**
```python
protocol = sim.SimProtocol(
    tree="((A:0.1,B:0.1):0.2,C:0.3);",
    root_seq_size=500,
    deletion_rate=0.03,
    insertion_rate=0.01,
    seed=424242
)
```

#### Methods

##### Sequence Configuration

```python
protocol.set_sequence_size(sequence_size: int) -> None
protocol.get_sequence_size() -> int
protocol.set_min_sequence_size(min_sequence_size: int) -> None
```

##### Random Seed

```python
protocol.set_seed(seed: int) -> None
protocol.get_seed() -> int
```

##### Insertion Rates

```python
# Uniform rate across all branches
protocol.set_insertion_rates(insertion_rate: float) -> None

# Per-branch rates
protocol.set_insertion_rates(insertion_rates: List[float]) -> None

# Query rates
protocol.get_insertion_rate(branch_num: int) -> float
protocol.get_all_insertion_rates() -> Dict[int, float]
```

##### Deletion Rates

```python
# Uniform rate across all branches
protocol.set_deletion_rates(deletion_rate: float) -> None

# Per-branch rates
protocol.set_deletion_rates(deletion_rates: List[float]) -> None

# Query rates
protocol.get_deletion_rate(branch_num: int) -> float
protocol.get_all_deletion_rates() -> Dict[int, float]
```

##### Insertion Length Distributions

```python
# Uniform distribution across all branches
protocol.set_insertion_length_distributions(insertion_dist: Distribution) -> None

# Per-branch distributions
protocol.set_insertion_length_distributions(insertion_dists: List[Distribution]) -> None

# Query distributions
protocol.get_insertion_length_distribution(branch_num: int) -> Distribution
protocol.get_all_insertion_length_distribution() -> Dict[int, Distribution]
```

##### Deletion Length Distributions

```python
# Uniform distribution across all branches
protocol.set_deletion_length_distributions(deletion_dist: Distribution) -> None

# Per-branch distributions
protocol.set_deletion_length_distributions(deletion_dists: List[Distribution]) -> None

# Query distributions
protocol.get_deletion_length_distribution(branch_num: int) -> Distribution
protocol.get_all_deletion_length_distribution() -> Dict[int, Distribution]
```

### Simulator

Main simulation engine that executes MSA generation.

#### Constructor

```python
Simulator(
    simProtocol: SimProtocol,
    simulation_type: SIMULATION_TYPE
)
```

**Parameters:**
- `simProtocol`: Configured `SimProtocol` object
- `simulation_type`: One of `SIMULATION_TYPE.{NOSUBS, DNA, PROTEIN}`

**Example:**
```python
simulator = sim.Simulator(
    protocol,
    simulation_type=sim.SIMULATION_TYPE.PROTEIN
)
```

#### Methods

##### Substitution Model Configuration

```python
simulator.set_replacement_model(
    model: MODEL_CODES,
    amino_model_file: pathlib.Path = None,
    model_parameters: List[float] = None,
    gamma_parameters_alpha: float = 1.0,
    gamma_parameters_categories: int = 1,
    invariant_sites_proportion: float = 0.0,
    site_rate_correlation: float = 0.0
) -> None
```

**Parameters:**
- `model`: Model code from `MODEL_CODES` enum (see [Substitution Models](#substitution-models))
- `amino_model_file`: Path to custom amino acid model file (for `MODEL_CODES.CUSTOM` only)
- `model_parameters`: Parameters for nucleotide models:
  - HKY: [freq_A, freq_C, freq_G, freq_T, kappa]
  - GTR: [freq_A, freq_C, freq_G, freq_T, rate_AC, rate_AG, rate_AT, rate_CG, rate_CT, rate_GT]
- `gamma_parameters_alpha`: Shape parameter for gamma rate heterogeneity (default: 1.0)
- `gamma_parameters_categories`: Number of discrete rate categories (default: 1)
- `invariant_sites_proportion`: Fraction of invariant sites (0.0-1.0, default: 0.0)
- `site_rate_correlation`: Autocorrelation coefficient for adjacent sites (0.0-1.0, default: 0.0)

**Example:**
```python
# Protein model with gamma rate heterogeneity
simulator.set_replacement_model(
    model=MODEL_CODES.WAG,
    gamma_parameters_alpha=1.0,
    gamma_parameters_categories=4
)

# HKY nucleotide model with custom parameters
simulator.set_replacement_model(
    model=MODEL_CODES.HKY,
    model_parameters=[0.25, 0.25, 0.25, 0.25, 2.0],  # Equal frequencies, kappa=2.0
    gamma_parameters_alpha=0.5,
    gamma_parameters_categories=4,
    invariant_sites_proportion=0.1
)
```

##### Sequence Saving Options

```python
simulator.save_root_sequence() -> None
simulator.save_all_nodes_sequences() -> None
simulator.get_sequences_to_save() -> List[bool]
```

By default, only leaf sequences are saved. Use these methods to also save ancestral sequences.

**Example:**
```python
# Save internal node sequences in addition to leaves
simulator.save_all_nodes_sequences()
```

##### Site Rate Management

```python
simulator.save_rates(is_save: bool) -> None
simulator.get_rates() -> List[float]
```

Enable rate tracking to retrieve per-site rates. Note: `get_rates()` empties the rate list after retrieval.

**Example:**
```python
simulator.save_rates(True)
msa = simulator()
rates = simulator.get_rates()
```

##### Simulation Execution

```python
# Standard mode (returns single MSA)
msa = simulator() -> Msa

# Multiple replicates (returns list of MSAs)
msas = simulator.simulate(times: int) -> List[Msa]

# Low-memory mode (writes directly to file)
simulator.simulate_low_memory(output_file_path: pathlib.Path) -> None
```

**Example:**
```python
# Single simulation
msa = simulator()

# 100 replicates
for i, msa in enumerate(simulator.simulate(times=100)):
    msa.write_msa(f"replicate_{i:04d}.fasta")

# Large simulation (low memory)
simulator.simulate_low_memory(pathlib.Path("large_output.fasta"))
```

##### Component-Level Simulation

Advanced methods for separate indel and substitution generation:

```python
blocktree = simulator.gen_indels() -> BlockTreePython
substitutions = simulator.gen_substitutions(length: int) -> sequenceContainer
```

### Distributions

Distribution classes define indel length probabilities.

#### Base Class

```python
class Distribution:
    """Abstract base class for indel length distributions."""
```

#### ZipfDistribution

Power-law distribution commonly used for indel lengths.

```python
ZipfDistribution(p: float, truncation: int = 150)
```

**Parameters:**
- `p`: Zipf exponent (typical values: 1.5-2.0; lower = longer indels)
- `truncation`: Maximum indel length

**Example:**
```python
# Moderate-length indels
dist = ZipfDistribution(1.7, 50)

# Longer indels
dist = ZipfDistribution(1.5, 100)
```

#### GeometricDistribution

Exponentially decreasing distribution.

```python
GeometricDistribution(p: float, truncation: int = 150)
```

**Parameters:**
- `p`: Geometric parameter (0 < p < 1; higher = shorter indels)
- `truncation`: Maximum indel length

**Example:**
```python
dist = GeometricDistribution(0.3, 50)
```

#### PoissonDistribution

Poisson-based distribution.

```python
PoissonDistribution(p: float, truncation: int = 150)
```

**Parameters:**
- `p`: Poisson lambda parameter
- `truncation`: Maximum indel length

**Example:**
```python
dist = PoissonDistribution(5.0, 50)
```

#### CustomDistribution

User-defined distribution from probability vector.

```python
CustomDistribution(dist: List[float])
```

**Parameters:**
- `dist`: List of probabilities (must sum to 1.0); index i is probability of length i+1

**Example:**
```python
# Custom distribution: 50% length-1, 30% length-2, 20% length-3
dist = CustomDistribution([0.5, 0.3, 0.2])
```

### Tree

Phylogenetic tree representation.

#### Constructor

```python
Tree(input_str: str)
```

**Parameters:**
- `input_str`: Newick format string or file path containing tree

**Example:**
```python
# From string
tree = sim.Tree("((A:0.1,B:0.1):0.2,C:0.3);")

# From file
tree = sim.Tree("path/to/tree.nwk")
```

#### Methods

```python
tree.get_num_nodes() -> int
tree.get_num_leaves() -> int
```

### Msa

Multiple sequence alignment result object.

#### Methods

##### Querying

```python
msa.get_length() -> int
msa.get_num_sequences() -> int
```

##### Output

```python
# Print to stdout
msa.print_msa() -> None
msa.print_indels() -> None

# Get as string
msa_string = msa.get_msa() -> str

# Write to file
msa.write_msa(file_path: str) -> None
```

**Example:**
```python
msa = simulator()

print(f"Alignment length: {msa.get_length()}")
print(f"Number of sequences: {msa.get_num_sequences()}")

# Write to FASTA
msa.write_msa("alignment.fasta")

# Get as string for processing
msa_str = msa.get_msa()
```

## Substitution Models

### Nucleotide Models

Available via `MODEL_CODES` enum:

- `NUCJC`: Jukes-Cantor (no parameters)
- `HKY`: Hasegawa-Kishino-Yano
  - Parameters: [freq_A, freq_C, freq_G, freq_T, kappa]
  - Example: `[0.25, 0.25, 0.25, 0.25, 2.0]`
- `GTR`: General Time Reversible
  - Parameters: [freq_A, freq_C, freq_G, freq_T, rate_AC, rate_AG, rate_AT, rate_CG, rate_CT, rate_GT]
- `TAMURA92`: Tamura 1992 model

### Protein Models

#### General Models
- `WAG`: Whelan and Goldman 2001
- `LG`: Le and Gascuel 2008
- `JONES`: Jones-Taylor-Thornton 1992 (JTT)
- `DAYHOFF`: Dayhoff et al. 1978
- `AAJC`: Amino acid Jukes-Cantor

#### Specialized Models
- `MTREV24`: Mitochondrial protein evolution
- `CPREV45`: Chloroplast protein evolution
- `HIVB`, `HIVW`: HIV between-host and within-host models

#### Structure-Based Models
- `EX_BURIED`: Buried sites
- `EX_EXPOSED`: Exposed sites
- `EHO_EXTENDED`: Extended structure
- `EHO_HELIX`: Helical structure
- `EHO_OTHER`: Other structures

#### Custom Models
- `CUSTOM`: User-provided amino acid model
  - Requires `amino_model_file` parameter with path to model file

### Usage Example

```python
from msasim.sailfish import MODEL_CODES

# Protein model
simulator.set_replacement_model(
    model=MODEL_CODES.WAG,
    gamma_parameters_alpha=1.0,
    gamma_parameters_categories=4
)

# DNA model with custom parameters
simulator.set_replacement_model(
    model=MODEL_CODES.HKY,
    model_parameters=[0.3, 0.2, 0.2, 0.3, 3.0],
    gamma_parameters_alpha=0.5,
    gamma_parameters_categories=4,
    invariant_sites_proportion=0.2
)
```

## Usage Patterns

### Pattern 1: Substitutions-Only Simulation

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
simulator.set_replacement_model(
    model=MODEL_CODES.LG,
    gamma_parameters_alpha=1.0,
    gamma_parameters_categories=4
)

msa = simulator()
msa.write_msa("substitutions_only.fasta")
```

### Pattern 2: Indels-Only Simulation

```python
from msasim import sailfish as sim
from msasim.sailfish import ZipfDistribution

protocol = sim.SimProtocol(
    tree="path/to/tree.nwk",
    root_seq_size=1000,
    deletion_rate=0.03,
    insertion_rate=0.01,
    deletion_dist=ZipfDistribution(1.7, 50),
    insertion_dist=ZipfDistribution(1.7, 50),
    seed=42
)

# NOSUBS type: no substitutions
simulator = sim.Simulator(protocol, simulation_type=sim.SIMULATION_TYPE.NOSUBS)
msa = simulator()
msa.write_msa("indels_only.fasta")
```

### Pattern 3: Full Model with Indels and Substitutions

```python
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES, ZipfDistribution

protocol = sim.SimProtocol(
    tree="path/to/tree.nwk",
    root_seq_size=500,
    deletion_rate=0.09,
    insertion_rate=0.03,
    deletion_dist=ZipfDistribution(1.7, 50),
    insertion_dist=ZipfDistribution(1.7, 50),
    seed=42
)

simulator = sim.Simulator(protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)
simulator.set_replacement_model(
    model=MODEL_CODES.WAG,
    gamma_parameters_alpha=1.0,
    gamma_parameters_categories=4
)

msa = simulator()
msa.write_msa("full_model.fasta")
```

### Pattern 4: Batch Simulations with Single Initialization

```python
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES

# Initialize once with base seed
protocol = sim.SimProtocol(tree="path/to/tree.nwk", root_seq_size=500, seed=42)
simulator = sim.Simulator(protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)
simulator.set_replacement_model(model=MODEL_CODES.JTT)

# Run multiple replicates
# Internal RNG advances automatically for reproducibility
for i in range(1000):
    msa = simulator()
    msa.write_msa(f"replicate_{i:05d}.fasta")
```

### Pattern 5: Per-Branch Parameter Configuration

```python
from msasim import sailfish as sim
from msasim.sailfish import ZipfDistribution

tree = sim.Tree("((A:0.1,B:0.2):0.3,C:0.4);")
protocol = sim.SimProtocol(tree=tree, root_seq_size=500, seed=42)

# Different rates per branch
num_branches = tree.get_num_nodes() - 1
insertion_rates = [0.01, 0.02, 0.015, 0.01, 0.03]
deletion_rates = [0.02, 0.03, 0.025, 0.02, 0.04]

protocol.set_insertion_rates(insertion_rates=insertion_rates)
protocol.set_deletion_rates(deletion_rates=deletion_rates)

# Different distributions per branch
distributions = [
    ZipfDistribution(1.5, 50),
    ZipfDistribution(1.7, 50),
    ZipfDistribution(1.6, 50),
    ZipfDistribution(1.8, 50),
    ZipfDistribution(1.5, 50),
]
protocol.set_insertion_length_distributions(insertion_dists=distributions)
protocol.set_deletion_length_distributions(deletion_dists=distributions)

simulator = sim.Simulator(protocol, simulation_type=sim.SIMULATION_TYPE.DNA)
msa = simulator()
```

### Pattern 6: Parameter Sweep

```python
from msasim import sailfish as sim
from msasim.sailfish import MODEL_CODES
import numpy as np

protocol = sim.SimProtocol(tree="path/to/tree.nwk", root_seq_size=500, seed=42)
simulator = sim.Simulator(protocol, simulation_type=sim.SIMULATION_TYPE.PROTEIN)

# Sweep gamma alpha parameter
for alpha in np.linspace(0.1, 5.0, 50):
    simulator.set_replacement_model(
        model=MODEL_CODES.WAG,
        gamma_parameters_alpha=alpha,
        gamma_parameters_categories=4
    )
    msa = simulator()
    msa.write_msa(f"alpha_{alpha:.2f}.fasta")
```

## Advanced Features

### Gamma Rate Heterogeneity

Models among-site rate variation using a discrete gamma distribution:

```python
simulator.set_replacement_model(
    model=MODEL_CODES.WAG,
    gamma_parameters_alpha=1.0,      # Shape parameter (lower = more variation)
    gamma_parameters_categories=4    # Number of rate categories (typically 4-8)
)
```

Lower `gamma_parameters_alpha` values indicate greater rate heterogeneity among sites.

### Invariant Sites

Specify a proportion of sites that never change:

```python
simulator.set_replacement_model(
    model=MODEL_CODES.LG,
    invariant_sites_proportion=0.2  # 20% of sites are invariant
)
```

### Site Rate Correlation

Model autocorrelation of rates along the sequence:

```python
simulator.set_replacement_model(
    model=MODEL_CODES.WAG,
    site_rate_correlation=0.5  # Moderate correlation between adjacent sites
)
```

### Ancestral Sequence Reconstruction

Save internal node sequences in addition to leaf sequences:

```python
# Save all internal nodes
simulator.save_all_nodes_sequences()
msa = simulator()

# Save only root
simulator.save_root_sequence()
msa = simulator()
```

### Site Rate Tracking

Retrieve per-site substitution rates:

```python
simulator.save_rates(True)
msa = simulator()
rates = simulator.get_rates()

# rates[i] contains the rate for site i
```

### Low-Memory Mode

For large simulations, write directly to disk without holding MSA in memory:

```python
import pathlib

simulator.simulate_low_memory(pathlib.Path("large_alignment.fasta"))
```

Use when:
- Tree has > 100K sequences
- Alignment length > 10K sites
- Available RAM is limited

Rough memory usage estimate: `(num_sequences * alignment_length) / 300,000` MB

## Performance Considerations

### Memory Management

For large datasets:

1. Use `simulate_low_memory()` instead of standard mode
2. Delete MSA objects immediately after use:
   ```python
   msa = simulator()
   msa.write_msa("output.fasta")
   del msa
   ```
3. Run garbage collection between batches:
   ```python
   import gc
   gc.collect()
   ```

### Computational Efficiency

1. **Single initialization**: Initialize simulator once, call multiple times rather than re-initializing
2. **Appropriate tree sizes**: Simulation time scales linearly with tree size
3. **Sequence length**: Longer sequences require more memory but simulation time scales linearly
4. **Rate categories**: More gamma categories increase computation time

### Reproducibility

1. **Set seed explicitly**: Always specify seed for reproducible results
2. **Single initialization pattern**: Initialize once with base seed, call repeatedly
   - Internal RNG advances deterministically across calls
3. **Avoid re-initialization**: Re-initializing with same seed resets RNG state

**Example:**
```python
# Good: reproducible per-replicate results
simulator = build_simulator(seed=42)
for i in range(100):
    msa = simulator()  # Each replicate has reproducible output

# Bad: all replicates identical
for i in range(100):
    simulator = build_simulator(seed=42)
    msa = simulator()  # Same output every iteration
```

### Typical Performance

Approximate simulation times on modern hardware:

- 1K sequences, 10K sites: < 1 second
- 10K sequences, 30K sites: ~10 seconds
- 100K sequences, 30K sites: ~2 minutes
- 1M sequences, 30K sites: ~20 minutes

Performance varies with:
- Indel rates (higher rates = longer runtime)
- Number of gamma categories (more categories = longer runtime)
- Tree topology (more balanced trees = faster)
