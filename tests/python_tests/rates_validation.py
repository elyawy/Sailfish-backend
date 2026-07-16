from msasim import ReplacementModelSpec, SiteRateModelSpec , Simulator, SimProtocol
from msasim.constants import ALPHABET_CODES, MODEL_CODES, SITE_RATE_MODELS


simulation_protocol = SimProtocol("((A:0.1,B:0.05):0.2,C:0.05);")
simulation_protocol.set_seed(42)
simulation_protocol.set_sequence_size(100)
simulation_protocol.set_max_insertion_length(10)
simulation_protocol.set_site_rate_model(SITE_RATE_MODELS.INDEL_AWARE)
# time.sleep(3)

site_rate_model_spec = SiteRateModelSpec(gamma_alpha=0.5,
                                        gamma_categories=8,
                                        site_rate_correlation=0.99,
                                        indel_awareness=SITE_RATE_MODELS.INDEL_AWARE)

rep_model = ReplacementModelSpec(MODEL_CODES.NUCJC, ALPHABET_CODES.DNA, site_rate_model=site_rate_model_spec)

simulator = Simulator(simulation_protocol, rep_model)

simulator.save_rates(True)

# msa = simulator()

for _ in range(4):
    msa = simulator()
    rates = simulator.get_rates()
    rate_categories = simulator.get_rate_categories()

    for i in range(len(rate_categories)):
        print(f"Site {i}: Rate category {rate_categories[i]}, Rate {rates[i]}")

    print()