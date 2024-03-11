import sys, pathlib, time
sys.path.insert(0,str(pathlib.Path(".").resolve()))
from Sailfish import simulator as sim
print(sim.__file__)


length_insertions = sim.ZipfDistribution(1.5, 100)
length_deletions = sim.GeometricDistribution(0.3, 50)


rate_insertion = 0.02
rate_deletion = 0.03

tree = sim.Tree("/home/elyalab/Data/Bacillus_1ZARG/RAxML_tree.tree")
print(tree.get_num_nodes())

simulation_protocol = sim.SimProtocol("/home/elyalab/Data/Bacillus_1ZARG/RAxML_tree.tree")
simulation_protocol.set_seed(10)

simulation_protocol.set_insertion_length_distributions(length_insertions)
simulation_protocol.set_deletion_length_distributions(length_deletions)
simulation_protocol.set_insertion_rates(rate_insertion)
simulation_protocol.set_deletion_rates(rate_deletion)
simulation_protocol.set_sequence_size(1000)
time.sleep(3)

simulator = sim.Simulator(simulation_protocol)
msas = simulator.simulate(1000)
