We provide additional usage examples here. The following examples demonstrate how to use the `msasim` package for different simple simulation scenarios.

## Substitutions only (no indels)

Set both rates indel rates to zero:

```python
from msasim import SimProtocol, Simulator, ReplacementModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES

protocol = SimProtocol(
    tree="(A:0.5,B:0.5);",
    root_seq_size=100,
    insertion_rate=0.0, # set to zero to disable insertions
    deletion_rate=0.0, # set to zero to disable deletions
    seed=42,
)

rep_model = ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA)
simulator = Simulator(protocol, replacement_model=rep_model)
msa = simulator()
```

## Indels only (no substitutions)

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

## Gamma rate heterogeneity

Set the `site_rate_model` argument in the `ReplacementModelSpec` to a `SiteRateModelSpec` instance:

```python
from msasim import SimProtocol, Simulator, ReplacementModelSpec, SiteRateModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES

protocol = SimProtocol(tree="(A:0.5,B:0.5);", root_seq_size=100, seed=42)

rate_model = SiteRateModelSpec(gamma_alpha=1.0, gamma_categories=4) # here we specify a gamma distribution with alpha=1.0 and 4 discrete categories
rep_model = ReplacementModelSpec(
    model=MODEL_CODES.WAG,
    alphabet=ALPHABET_CODES.PROTEIN,
    site_rate_model=rate_model,
)
simulator = Simulator(protocol, replacement_model=rep_model)
msa = simulator()
```

## Batch simulations

Simulate multiple replicates in a loop. The internal random number generator (RNG) advances deterministically across calls, so each replicate will be different but reproducible if the same seed is used. Single evolutionary model is used to simulate all replicates.

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

## Low-memory mode for large simulations

During simulation, the entire multiple sequence alignment (MSA) is not stored in memory. Instead, each replicate is written to disk immediately during simulation. This allows for simulating very large alignments (1M+ sites) without running out of memory.

```python
import pathlib
from msasim import SimProtocol, Simulator, ReplacementModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES

protocol = SimProtocol(tree="large_tree.nwk", root_seq_size=10000, seed=42)
rep_model = ReplacementModelSpec(model=MODEL_CODES.NUCJC, alphabet=ALPHABET_CODES.DNA)
simulator = Simulator(protocol, replacement_model=rep_model)

simulator.simulate(output_path=pathlib.Path("large_alignment.fasta")) # → writes large_alignment_replicate_1.fasta
```
