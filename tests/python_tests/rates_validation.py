from msasim import sailfish as sim



simulation_protocol = sim.SimProtocol("((A:0.1,B:0.05):0.02,C:0.05);")
simulation_protocol.set_seed(42)
simulation_protocol.set_sequence_size(100)
simulation_protocol.set_max_insertion_length(10)
simulation_protocol.set_site_rate_model(sim.SITE_RATE_MODELS.INDEL_AWARE)
# time.sleep(3)

simulator = sim.Simulator(simulation_protocol, simulation_type=sim.SIMULATION_TYPE.DNA)
simulator.set_replacement_model(model=sim.MODEL_CODES.NUCJC,
                                gamma_parameters_alpha=1.0,
                                gamma_parameters_categories=8,
                                site_rate_correlation=0.5)

simulator.save_rates(True)

# msa = simulator()

for _ in range(4):
    msa = simulator()
    rates = simulator.get_rates()
    rate_categories = simulator.get_rate_categories()

    for i in range(len(rate_categories)):
        print(f"Site {i}: Rate category {rate_categories[i]}, Rate {rates[i]}")

    print()