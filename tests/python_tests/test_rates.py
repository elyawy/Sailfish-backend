import sys, pathlib, time
sys.path.insert(0,str(pathlib.Path(".").resolve()))
from msasim import sailfish as sim
print(sim.__file__)



simulation_protocol = sim.SimProtocol("((A:0.1,B:0.05):0.02,C:0.05);")
simulation_protocol.set_seed(1)
simulation_protocol.set_sequence_size(10)
# time.sleep(3)

simulator = sim.Simulator(simulation_protocol, simulation_type=sim.SIMULATION_TYPE.DNA)
simulator.set_replacement_model(model=sim.MODEL_CODES.NUCJC,
                                gamma_parameters_alpha=2.0,
                                gamma_parameters_catergories=50)

simulator.save_rates(True)

msa = simulator()
print(simulator.get_rates()) # after calling this function the list is emptied
print(simulator.get_rates()) # this is empty []

msa = simulator() # filled rates back with new simulation.
msa.print_msa()
print(simulator.get_rates())