import sys, pathlib
from _Sailfish import Tree



tree = Tree(f"tests/trees/normalbranches_nLeaves10.treefile")

root_node = tree.root


assert root_node.distance_to_father() == -1.0
sum_nodes = 0
node_list = [root_node]
node_distances = []
# print(node_list)

for node in node_list:
    for son in node.sons:
        node_list.append(son)
        node_distances.append(son.distance_to_father())




print(sum(node_distances))