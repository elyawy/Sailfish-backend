import pickle
from _Sailfish import Tree, SimProtocol, DiscreteDistribution


tree = Tree("tests/trees/normalbranches_nLeaves10.treefile", True)


with open("/home/elyalab/Data/simulator_tests/tree.pickle", 'wb') as f:
    pickle.dump(tree, f)

with open("/home/elyalab/Data/simulator_tests/tree.pickle", 'rb') as f:
    pickled_tree = pickle.load(f)

print(tree.root.sons[1].distance_to_father(), pickled_tree.root.sons[1].distance_to_father())
print(tree.root.sons[1].distance_to_father(), pickled_tree.root.sons[1].distance_to_father())

simprot = SimProtocol(tree)
simprot.set_seed(10)
simprot.set_insertion_rates([0.02]* (tree.num_nodes - 1))
simprot.set_deletion_rates([0.04]* (tree.num_nodes -1))
dist = DiscreteDistribution([0.5,0.5])
simprot.set_deletion_length_distributions([dist]*(tree.num_nodes - 1))
simprot.set_insertion_length_distributions([dist]*(tree.num_nodes - 1))



with open("/home/elyalab/Data/simulator_tests/simprot.pickle", 'wb') as f:
    pickle.dump(simprot, f)

with open("/home/elyalab/Data/simulator_tests/simprot.pickle", 'rb') as f:
    pickled_simprot = pickle.load(f)

print(pickled_simprot.get_seed())
print(pickled_simprot.get_deletion_rate(5))
print(pickled_simprot.get_insertion_length_distribution(3))

# print(simprot.get_deletion_length_distribution(0), pickled_simprot.get_deletion_length_distribution(0))
