"""
Recipe: Continuous Gamma Rate Heterogeneity
============================================
Sailfish supports up to 254 rate categories per simulation. For sequences
longer than 254 sites, we split the simulation into chunks of up to 254 sites,
each with rates sampled fresh from a continuous gamma distribution, then
concatenate the resulting MSAs.

Requirements: numpy
"""

import math
import pathlib
import numpy as np

from msasim import SimProtocol, Simulator, ReplacementModelSpec, SiteRateModelSpec
from msasim import MODEL_CODES, ALPHABET_CODES

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

TREE = "(A:0.1,(B:0.05,C:0.05):0.05);"
ROOT_SEQ_SIZE = 1000       # total alignment length
GAMMA_ALPHA = 0.5          # shape parameter of the continuous gamma
MODEL = MODEL_CODES.JONES
ALPHABET = ALPHABET_CODES.PROTEIN
SEED = 1
CHUNK_SIZE = 254           # max rate categories per simulation

OUTPUT = pathlib.Path("continuous_gamma_output.fasta")


# ---------------------------------------------------------------------------
# Setup — build simulator once, reuse across chunks
# ---------------------------------------------------------------------------
 
rng = np.random.default_rng(SEED)
 
# Placeholder rate model — will be updated per chunk
initial_rates = rng.gamma(shape=GAMMA_ALPHA, scale=1.0 / GAMMA_ALPHA, size=CHUNK_SIZE)
rate_model = SiteRateModelSpec(
    free_rates=initial_rates.tolist(),
    free_rate_weights=[1.0 / CHUNK_SIZE] * CHUNK_SIZE,
)
rep_model = ReplacementModelSpec(model=MODEL, alphabet=ALPHABET, site_rate_model=rate_model)
 
protocol = SimProtocol(
    tree=TREE,
    root_seq_size=CHUNK_SIZE,
    insertion_rate=0.0,
    deletion_rate=0.0,
    seed=SEED,
)
protocol.set_cache_branch_probs(True)
 
simulator = Simulator(protocol, replacement_model=rep_model)
 
current_categories_size = None
 
# ---------------------------------------------------------------------------
# Simulate chunks
# ---------------------------------------------------------------------------
 
n_chunks = math.ceil(ROOT_SEQ_SIZE / CHUNK_SIZE)
chunk_sizes = [CHUNK_SIZE] * (n_chunks - 1) + [ROOT_SEQ_SIZE - CHUNK_SIZE * (n_chunks - 1)]
 
msas = []
for chunk_size in chunk_sizes:
    # Sample one rate per site from the continuous gamma
    rates = rng.gamma(shape=GAMMA_ALPHA, scale=1.0 / GAMMA_ALPHA, size=chunk_size)
 
    # Update the free rates for this chunk (cheap — does not rebuild the model)
    chunk_rate_model = SiteRateModelSpec(
        free_rates=rates.tolist(),
        free_rate_weights=[1.0 / chunk_size] * chunk_size,
    )
    simulator.set_replacement_model(
        ReplacementModelSpec(model=MODEL, alphabet=ALPHABET, site_rate_model=chunk_rate_model)
    )
 
    # Update sequence size and categories only when chunk size changes (i.e. the last chunk)
    if chunk_size != current_categories_size:
        protocol.set_sequence_size(chunk_size)
        simulator.set_per_site_rate_categories(list(range(chunk_size)))
        current_categories_size = chunk_size
 
    msas.append(simulator())
 
# ---------------------------------------------------------------------------
# Merge and write
# ---------------------------------------------------------------------------
 
partition_dicts = []
for msa in msas:
    rows = {}
    for i in range(msa.get_num_sequences()):
        name, seq = msa.get_msa_row(i).split("\n", 1)
        rows[name.lstrip(">")] = seq
    partition_dicts.append((rows, msa.get_length()))
 
taxa = sorted(partition_dicts[0][0].keys())
with open(OUTPUT, "w") as f:
    for taxon in taxa:
        full_seq = "".join(rows.get(taxon, "-" * length) for rows, length in partition_dicts)
        f.write(f">{taxon}\n{full_seq}\n")
 
print(f"Done — wrote {ROOT_SEQ_SIZE}-site MSA to {OUTPUT}")
