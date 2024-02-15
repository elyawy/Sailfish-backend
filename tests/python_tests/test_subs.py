import pathlib, time
from itertools import repeat
import numpy as np


from _Sailfish import Tree, Simulator, SimProtocol, Msa, DiscreteDistribution, modelFactory, alphabetCode, modelCode


# tree = Tree(f"{pathlib.Path.home()}/Data/yeast/RAxML_tree.tree")
tree = Tree(f"tests/trees/normalbranches_nLeaves10.treefile")

root_node = tree.root

protocol = SimProtocol(tree)


rates_i = list(repeat(0.01, tree.num_nodes))
rates_d = list(repeat(0.01, tree.num_nodes))

a_param = 1.08

truncation = 50 + 1

area_under = sum(i**(-a_param) for i in range(1, truncation))
# print(area_under)


len_dist = [(i**(-a_param))/area_under for i in range(1,truncation)]
# len_dist = [(i**(-a_param)) for i in range(1,51)]
# print(len_dist)
# len_dist.append(((50**(-a_param)) + tail)/area_under)
# print(len_dist)
# print(sum(len_dist))
# print(sum(i*len_dist[i-1] for i in range(1,51)))

# exit(1)
dist = DiscreteDistribution(len_dist)
rand_seed = int(str(time.time_ns())[-8:])
print(rand_seed)
dist.set_seed(rand_seed)


# freq_table = {i:0 for i in range(1,51)}
# for i in range(100000):
#     freq_table[dist.draw_sample()] += 1

# sum_dist = 0
# for i in range(1,51):
#     sum_dist += freq_table[i]/100000
#     print(freq_table[i]/100000)
# print(sum_dist)

all_dists = list(repeat(dist, tree.num_nodes))

# print(all_dists)


protocol.set_sequence_size(100)
protocol.set_insertion_rates(rates_i)
protocol.set_deletion_rates(rates_d)
protocol.set_insertion_length_distributions(all_dists)
protocol.set_deletion_length_distributions(all_dists)

protocol.set_seed(rand_seed)

sim = Simulator(protocol)


mFac = modelFactory(tree)

mFac.set_alphabet(alphabetCode.AMINOACID)
mFac.set_replacement_model(modelCode.WAG)
mFac.set_gamma_parameters(0.5, 4)


min_lengths = []
max_lengths = []

msas = []
substitution_list = []


for i in range(10):
    blockmap = sim.gen_indels()
    msa = Msa(blockmap, root_node)
    substitutions = sim.gen_substitutions(mFac, msa.length())
    substitution_list.append(substitutions)
    msa.fill_substitutions(substitutions)
    msas.append(msa)

# Hack to catch cout output as string
# import os, io, sys, tempfile

# def redirect_stdout(to_fd, stdout):
#     """Redirect stdout to the given file descriptor."""
#     sys.stdout.close()
#     os.dup2(to_fd, stdout)
#     sys.stdout = io.TextIOWrapper(os.fdopen(stdout, 'wb'))


# stdout = sys.stdout.fileno()
# saved = os.dup(stdout)
# with tempfile.TemporaryFile() as f:
# # redirect stdout to a buffer
#     redirect_stdout( f.fileno(), stdout )
#     msas[4].print_msa()
#     redirect_stdout( saved, stdout )
#     f.seek(0)
#     fasta = f.read().decode().split("\n")

# names = [f">T{i}\n" for i in range(1,11)]
# fasta = list(zip(names, fasta))
# fasta = [i + j for i,j in fasta]
# fasta = "\n".join(fasta)
# print(fasta)