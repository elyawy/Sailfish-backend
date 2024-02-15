import pathlib, time
from itertools import repeat
import numpy as np


from _Sailfish import Tree, Simulator, SimProtocol, Msa, DiscreteDistribution


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


protocol.set_sequence_size(1000)
protocol.set_insertion_rates(rates_i)
protocol.set_deletion_rates(rates_d)
protocol.set_insertion_length_distributions(all_dists)
protocol.set_deletion_length_distributions(all_dists)

protocol.set_seed(rand_seed)

sim = Simulator(protocol)

min_lengths = []
max_lengths = []

for i in range(10000):
    blockmap = sim.gen_indels()
    msa = Msa(blockmap, root_node)
    seq_lengths = []
    for seq in msa.get_msa():
        current_arr = np.array(seq)
        # print("".join([item*"A" if item > 0 else (-item)*"-" for item in seq]))
        seq_lengths.append(current_arr[current_arr > 0].sum())
    # print()
    min_lengths.append(min(seq_lengths))
    max_lengths.append(max(seq_lengths))

print("min mean:", sum(min_lengths)/len(min_lengths))
print("max mean:", sum(max_lengths)/len(max_lengths))

print(max(max_lengths))

from matplotlib import pyplot as plt

plt.hist(min_lengths, bins=50)
plt.show()

plt.hist(max_lengths, bins=50)
plt.show()
