import sys, pathlib, time
sys.path.insert(0,str(pathlib.Path(".").resolve()))
import psutil
from Sailfish import simulator as sim
print(sim.__file__)


length_insertions = sim.ZipfDistribution(1.7, 50)
length_deletions = sim.ZipfDistribution(1.7, 50)

rate_insertion = 0.03
rate_deletion = 0.09

# tree = sim.Tree("/home/elyalab/Data/Bacillus_1ZARG/RAxML_tree.tree")
# print(tree.get_num_nodes())

trees_path = pathlib.Path("tests/trees").resolve()

trees_map = {
    "10": trees_path / "normalbranches_nLeaves10.treefile",
    "100": trees_path / "normalbranches_nLeaves100.treefile",
    "1k": trees_path / "normalbranches_nLeaves1000.treefile",
    "5k": trees_path / "normalbranches_nLeaves5000.treefile",
    # "10k": trees_path / "normalbranches_nLeaves10000.treefile"
    # "100k": trees_path / "normalbranches_nLeaves100000.treefile"
}

def init_protocol(number_of_species) -> sim.Simulator:
    simulation_protocol = sim.SimProtocol(str(trees_map[number_of_species]))
    simulation_protocol.set_seed(100)

    simulation_protocol.set_insertion_length_distributions(length_insertions)
    simulation_protocol.set_deletion_length_distributions(length_deletions)
    simulation_protocol.set_insertion_rates(rate_insertion)
    simulation_protocol.set_deletion_rates(rate_deletion)
    simulation_protocol.set_sequence_size(30000)
    # time.sleep(3)

    simulator = sim.Simulator(simulation_protocol, simulation_type=sim.SIMULATION_TYPE.DNA)
    simulator.set_replacement_model(model=sim.MODEL_CODES.NUCJC)

    return simulator

def time_me(func):
    def wrapper(*args, **kwargs):
        tic = time.perf_counter()
        res = func(*args, **kwargs)
        toc = time.perf_counter()
        print(f"function {func.__name__} took {toc - tic:0.10f} seconds with {num_sequences} sequences")
        return res

    return wrapper



for num_sequences in trees_map.keys():
    simulator = init_protocol(num_sequences)
    msa = time_me(simulator.simulate)()
    process = psutil.Process()
    print(process.memory_info().rss / 1024**3)  # in bytes 
    # msa.write_msa(f"simulator_tests/msa_{num_sequences}.fasta") 
