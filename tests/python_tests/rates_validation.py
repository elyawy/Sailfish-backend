from msasim import ReplacementModelSpec, SiteRateModelSpec , Simulator, SimProtocol
from msasim.constants import ALPHABET_CODES, MODEL_CODES, SITE_RATE_MODELS
from msasim.distributions import ZipfDistribution

simulation_protocol = SimProtocol("tests/trees/normalbranches_nLeaves100.treefile")
simulation_protocol.set_seed(42)
simulation_protocol.set_sequence_size(1000)
simulation_protocol.set_insertion_rates(1.0)
simulation_protocol.set_deletion_rates(0.0)
simulation_protocol.set_insertion_length_distributions(ZipfDistribution(1.7, 50))
simulation_protocol.set_deletion_length_distributions(ZipfDistribution(1.7, 50))
# time.sleep(3)

# site_rate_model_spec = SiteRateModelSpec(gamma_alpha=0.5,
#                                         gamma_categories=8,
#                                         # site_rate_correlation=0.9,
#                                         indel_awareness=SITE_RATE_MODELS.SIMPLE)
# rep_model = ReplacementModelSpec(MODEL_CODES.NUCJC, ALPHABET_CODES.DNA, site_rate_model=site_rate_model_spec)
simulation_protocol._set_site_rate_model(SITE_RATE_MODELS.INDEL_AWARE)
simulator = Simulator(simulation_protocol)

# simulator.save_rates(True)
msa = simulator()

# for _ in range(4):
#     msa = simulator()
#     rates = simulator.get_rates()
#     rate_categories = simulator.get_rate_categories()

#     # print the Site numbers horizontally with the corresponding rate category underneath and rate underneath the category
#     # aligned correctly with the site number and category, rate have 4 significant figures
#     print("Site: ", end="")
#     for site in range(len(rates)):
#         print(f"{site:>6}", end="")

#     print("\nCat:  ", end="")
#     for cat in rate_categories:
#         print(f"{cat:>6}", end="")
    
#     print("\nRate: ", end="")
#     for rate in rates:
#         print(f"{rate:>6.2f}", end="")

#     print("\n")