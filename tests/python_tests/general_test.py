
# import sys
# sys.path.append("/home/elyalab/Dev/projects/Sailfish-backend/Sailfish/")
from Sailfish import simulator as sim
import numpy as np

TRUNCATION_NUM = 20
NUM_OF_SAMPLES = 10000

geodist = sim.GeometricDistribution(0.3, TRUNCATION_NUM)

samples = [geodist.draw_sample() for i in range(NUM_OF_SAMPLES)]

frq = {i:0 for i in range(1, TRUNCATION_NUM + 1)}
for x in samples:
    frq[x] = frq.get(x, 0) + 1

print(f'geometric dist: ')
print(frq)

possion = sim.PoissonDistribution(10, TRUNCATION_NUM)

samples = [possion.draw_sample() for i in range(NUM_OF_SAMPLES)]

frq = {i:0 for i in range(1, TRUNCATION_NUM + 1)}
for x in samples:
    frq[x] = frq.get(x, 0) + 1

print(f'possion dist: ')
print(frq)


tree = sim.Tree("(A:0.01,B:0.5,C:0.03);")

print(tree.get_num_nodes())


tree = sim.Tree("/home/elyalab/Dev/projects/Sailfish-backend/tests/trees/normalbranches_nLeaves100.treefile")

print(tree.get_num_nodes())

protocol = sim.SimProtocol(tree=tree)

protocol.set_insertion_length_distributions(geodist)
protocol.set_deletion_length_distributions(possion)

protocol.set_insertion_rates(0.1)
protocol.set_deletion_rates(0.1)
protocol.set_sequence_size(100)

print(protocol.get_all_deletion_length_distribution())
print(protocol.get_all_insertion_length_distribution())
print(protocol.get_all_deletion_rates())
print(protocol.get_all_insertion_rates())

protocol.set_seed(2)

simantov = sim._Sailfish.Simulator(protocol.sim)

blocktree = simantov.gen_indels()

msa = sim._Sailfish.Msa(blocktree, tree._tree.root)

msa.print_indels()