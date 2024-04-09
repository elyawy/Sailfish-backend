import sys, pathlib, time
sys.path.insert(0,str(pathlib.Path(".").resolve()))

from Sailfish import simulator

protocol =  simulator.SimProtocol(tree="tests/trees/normalbranches_nLeaves1000.treefile", 
                                  deletion_rate=0.01,
                                  insertion_rate=0.01, seed=42, root_seq_size=30000)


sim = simulator.Simulator(protocol, simulation_type=simulator.SIMULATION_TYPE.DNA)

sim()