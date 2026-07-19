"""
Recipe: Hypermutability
=========================
This recipe demonstrates how to simulate a sequence alignment 
with a mixture of standard and hypermutable substitution models.
The model here is designed to mimic viral hypermutability, where certain sites
are prone to rapid mutation (e.g., G -> T transitions).
The simulation uses a mixture of two models:
1. A standard GTR model for the majority of sites.

"""

import pathlib
from msasim import (
    SimProtocol, 
    ReplacementModelSpec, 
    MODEL_CODES, 
    ALPHABET_CODES
)
from msasim.advanced import MixtureModel, Mixture

# ---------------------------------------------------------
# 1. Establish the Simulation Protocol
# ---------------------------------------------------------
# Set up your tree and ensure insertion/deletion rates are 0.0 
# if you want a substitutions-only model.
protocol = SimProtocol(
    tree="(((A:0.1, B:0.1):0.05, C:0.15):0.05, D:0.2);",
    root_seq_size=500,  # Total length of the simulated alignment
    insertion_rate=0.0,
    deletion_rate=0.0,
    seed=333
)


# ---------------------------------------------------------
# 2. Define the Standard Background Model (e.g., GTR)
# ---------------------------------------------------------
# Standard GTR parameters: [freqA, freqC, freqG, freqT, rAC, rAG, rAT, rCG, rCT, rGT]
background_gtr = ReplacementModelSpec(
    model=MODEL_CODES.GTR,
    alphabet=ALPHABET_CODES.DNA,
    model_parameters=[
        0.25, 0.25, 0.25, 0.25,  # Equal state frequencies
        1.0,  3.0,  1.0,  1.0,  5.0,  1.0  # Normalized rate ratios
    ]
)

# ---------------------------------------------------------
# 3. Define the Non-Reversible Hypermutable Model
# ---------------------------------------------------------
# For a DNA alphabet (size 4), the NONREVERSIBLE model requires a 20-parameter list:
# - A 4x4 instantaneous Q-rate matrix (row-by-row: A, C, G, T)
# - Followed by the 4 equilibrium frequencies.
# We will construct a matrix where G -> T is scaled 100-fold.

Q_matrix = [
    # Row A: (A->A is auto-calculated, A->C, A->G, A->T)
    -0.4,   0.1,   0.2,   0.1, 
    # Row C: (C->A, C->C, C->G, C->T)
    0.1,  -0.4,   0.1,   0.2,
    # Row G: (G->A, G->C, G->G, G->T is HYPERMUTABLE)
    0.1,  0.1,  -20.2,  20.0, # Massive boost to G -> T
    # Row T: (T->A, T->C, T->G, T->T)
    0.1,  0.2,  0.1,  -0.4,  
]
eq_freqs = [0.25, 0.25, 0.25, 0.25]

hypermutable_model = ReplacementModelSpec(
    model=MODEL_CODES.NONREVERSIBLE,
    alphabet=ALPHABET_CODES.DNA,
    model_parameters=Q_matrix + eq_freqs
)

# ---------------------------------------------------------
# 4. Compose the Mixture
# ---------------------------------------------------------
# We assign 98% of the sites to the standard background model,
# and 2% of the sites to be stochastically hypermutable.
models = [background_gtr, hypermutable_model]
weights = [0.75, 0.25] 

mixture_setup = MixtureModel(
    name="viral_hypermutability_mixture",
    replacement_models=models,
    model_weights=weights,
    indel_model=protocol
)

# ---------------------------------------------------------
# 5. Simulate the Alignment
# ---------------------------------------------------------
simulator = Mixture(mixture_setup)
protocol.get_sim_context().set_save_root()

msa = simulator.simulate()

# Output the simulated sequences to a file
msa.write_msa("hypermutable_sim.fasta")
print("Simulation complete! Alignment length:", msa.get_length())