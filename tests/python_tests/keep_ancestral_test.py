
from msasim import MODEL_CODES, ZipfDistribution, SimProtocol, Simulator
from msasim.substitutions import ReplacementModelSpec, SiteRateModelSpec

ROOT_SEQUENCE_LENGTH = 150

sim_protocol = SimProtocol("((A:0.1,B:1.0):0.05,(C:0.2):0.01);",
                               deletion_rate=0.01,
                               insertion_rate=0.01,
                               deletion_dist=ZipfDistribution(1.08, 50),
                               insertion_dist=ZipfDistribution(1.08, 50),
                               seed=9)

sim_protocol.set_sequence_size(ROOT_SEQUENCE_LENGTH)

site_rate_model = SiteRateModelSpec(gamma_alpha=1.0,
                                   gamma_categories=4)

replacement_model_spec = ReplacementModelSpec(model=MODEL_CODES.WAG, 
                                              site_rate_model=site_rate_model)

simulation = Simulator(sim_protocol, replacement_model_spec)

simulation.save_root_sequence()

msa = simulation()
msa.print_msa()
