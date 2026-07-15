#include <emscripten/bind.h>
#include <memory>
#include <vector>
#include <string>

#include "../libs/pcg/pcg_random.hpp"
// #include "../libs/splitmix/Splitmix64.h"
#include "../libs/Phylolib/includes/gammaDistribution.h"
#include "./IndelSimulator.h"
#include "./SubstitutionSimulator.h"

using namespace emscripten;
using SelectedRNG = pcg64_fast;

EMSCRIPTEN_BINDINGS(_Sailfish) {

    // -------------------------------------------------------------------------
    // Opaque vector types
    // value_type extracts the element type from the typedef, no need to hardcode
    // -------------------------------------------------------------------------
    register_vector<SparseSequenceContainer::value_type>("SparseSequenceContainer");
    register_vector<size_t>("IntVector");
    register_vector<SparseMSA::value_type>("SparseMSA");
    register_vector<double>("DoubleVector");
    register_vector<std::vector<double>>("DoubleMatrix");
    register_vector<Event>("EventSequence");
    register_vector<EventSequence>("EventMap");
    register_vector<tree::TreeNode*>("NodeVector");
    // -------------------------------------------------------------------------
    // Tree
    // -------------------------------------------------------------------------
    class_<tree::TreeNode>("node")
        .property("sons",       &tree::TreeNode::getSons, allow_raw_pointers())
        .property("num_leaves", &tree::TreeNode::getNumberLeaves)
        .property("name",       &tree::TreeNode::name)
        .function("distance_to_father", &tree::TreeNode::dis2father);

    class_<tree>("Tree")
        // Constructor accepts newick string OR file path — same as C++ side
        .constructor<const std::string&, bool>()
        .property("num_nodes", &tree::getNodesNum)
        .property("root",      &tree::getRoot, allow_raw_pointers());

    // -------------------------------------------------------------------------
    // Enums
    // -------------------------------------------------------------------------
    enum_<SiteRateModel>("SiteRateModel")
        .value("SIMPLE",      SiteRateModel::SIMPLE)
        .value("INDEL_AWARE", SiteRateModel::INDEL_AWARE);

    enum_<modelCode>("modelCode")
        .value("NUCJC",          modelCode::NUCJC)
        .value("AAJC",           modelCode::AAJC)
        .value("GTR",            modelCode::GTR)
        .value("HKY",            modelCode::HKY)
        .value("TAMURA92",       modelCode::TAMURA92)
        .value("CPREV45",        modelCode::CPREV45)
        .value("DAYHOFF",        modelCode::DAYHOFF)
        .value("JONES",          modelCode::JONES)
        .value("MTREV24",        modelCode::MTREV24)
        .value("WAG",            modelCode::WAG)
        .value("HIVB",           modelCode::HIVB)
        .value("HIVW",           modelCode::HIVW)
        .value("LG",             modelCode::LG)
        .value("EMPIRICODON",    modelCode::EMPIRICODON)
        .value("EX_BURIED",      modelCode::EX_BURIED)
        .value("EX_EXPOSED",     modelCode::EX_EXPOSED)
        .value("EHO_EXTENDED",   modelCode::EHO_EXTENDED)
        .value("EHO_HELIX",      modelCode::EHO_HELIX)
        .value("EHO_OTHER",      modelCode::EHO_OTHER)
        .value("EX_EHO_BUR_EXT", modelCode::EX_EHO_BUR_EXT)
        .value("EX_EHO_BUR_HEL", modelCode::EX_EHO_BUR_HEL)
        .value("EX_EHO_BUR_OTH", modelCode::EX_EHO_BUR_OTH)
        .value("EX_EHO_EXP_EXT", modelCode::EX_EHO_EXP_EXT)
        .value("EX_EHO_EXP_HEL", modelCode::EX_EHO_EXP_HEL)
        .value("EX_EHO_EXP_OTH", modelCode::EX_EHO_EXP_OTH)
        .value("CUSTOM",         modelCode::CUSTOM);

    enum_<event>("IndelEventType")
        .value("INSERTION", event::INSERTION)
        .value("DELETION",  event::DELETION);

    // -------------------------------------------------------------------------
    // DiscreteDistribution
    // -------------------------------------------------------------------------
    class_<DiscreteDistribution>("DiscreteDistribution")
        .constructor<std::vector<double>>();

    // -------------------------------------------------------------------------
    // GammaDistribution
    // -------------------------------------------------------------------------
    class_<gammaDistribution>("GammaDistribution")
        .constructor<MDOUBLE, int>()
        .function("getAllRates", optional_override([](const gammaDistribution& g) {
            std::vector<MDOUBLE> result;
            for (size_t i = 0; i < g.categories(); ++i)
                result.push_back(g.rates(i));
            return result;
        }))
        .function("getAllRatesProb", optional_override([](const gammaDistribution& g) {
            std::vector<MDOUBLE> result;
            for (size_t i = 0; i < g.categories(); ++i)
                result.push_back(g.ratesProb(i));
            return result;
        }));

    // -------------------------------------------------------------------------
    // SimulationProtocol
    // -------------------------------------------------------------------------
    class_<SimulationProtocol>("SimProtocol")
        .constructor<size_t>()
        .function("set_sequence_size",                   &SimulationProtocol::setSequenceSize)
        .function("get_sequence_size",                   &SimulationProtocol::getSequenceSize)
        .function("set_insertion_rates",                 &SimulationProtocol::setInsertionRates)
        .function("get_insertion_rate",                  &SimulationProtocol::getInsertionRate)
        .function("set_deletion_rates",                  &SimulationProtocol::setDeletionRates)
        .function("get_deletion_rate",                   &SimulationProtocol::getDeletionRate)
        // Accepts a single DiscreteDistribution applied uniformly to all edges.
        // (Exposing vector<DiscreteDistribution*> directly is not feasible from JS.)
        .function("set_insertion_length_distribution_uniform", optional_override([](SimulationProtocol& self, DiscreteDistribution& dist, size_t numEdges) {
            std::vector<DiscreteDistribution*> dists(numEdges, &dist);
            self.setInsertionLengthDistributions(dists);
        }))
        .function("get_insertion_length_distribution",   &SimulationProtocol::getInsertionDistribution, allow_raw_pointers())
        .function("set_deletion_length_distribution_uniform", optional_override([](SimulationProtocol& self, DiscreteDistribution& dist, size_t numEdges) {
            std::vector<DiscreteDistribution*> dists(numEdges, &dist);
            self.setDeletionLengthDistributions(dists);
        }))
        .function("get_deletion_length_distribution",    &SimulationProtocol::getDeletionDistribution, allow_raw_pointers())
        .function("set_insertion_length_distributions", optional_override([](SimulationProtocol& self, emscripten::val distArray, size_t numEdges) {
            std::vector<DiscreteDistribution*> dists;
            dists.reserve(numEdges);
            for (size_t i = 0; i < numEdges; ++i)
                dists.push_back(distArray[i].as<DiscreteDistribution*>(allow_raw_pointers()));
            self.setInsertionLengthDistributions(dists);
        }))
        .function("set_deletion_length_distributions", optional_override([](SimulationProtocol& self, emscripten::val distArray, size_t numEdges) {
            std::vector<DiscreteDistribution*> dists;
            dists.reserve(numEdges);
            for (size_t i = 0; i < numEdges; ++i)
                dists.push_back(distArray[i].as<DiscreteDistribution*>(allow_raw_pointers()));
            self.setDeletionLengthDistributions(dists);
        }))
        .function("set_minimum_sequence_size",           &SimulationProtocol::setMinSequenceSize)
        .function("get_minimum_sequence_size",           &SimulationProtocol::getMinSequenceSize)
        .function("set_site_rate_model",                 &SimulationProtocol::setSiteRateModel)
        .function("get_site_rate_model",                 &SimulationProtocol::getSiteRateModel)
        .function("set_max_insertion_length",            &SimulationProtocol::setMaxInsertionLength)
        .function("get_max_insertion_length",            &SimulationProtocol::getMaxInsertionLength);

    // -------------------------------------------------------------------------
    // CategorySampler (opaque handle — only passed around, not constructed)
    // -------------------------------------------------------------------------
    class_<CategorySampler>("CategorySampler");

    // -------------------------------------------------------------------------
    // modelFactory
    //
    // set_site_rate_model: transition_matrix was optional in pybind11.
    // Embind has no default args, so expose two functions:
    //   set_site_rate_model(rates, probs, matrix)
    //   set_site_rate_model_no_matrix(rates, probs)          <- pass [] from JS
    //
    // get_rate_category_sampler: max_path_length defaulted to 0.
    //   get_rate_category_sampler(max_path_length)
    //   get_rate_category_sampler_default()                  <- uses 0
    // -------------------------------------------------------------------------
    class_<modelFactory>("modelFactory")
        .constructor<>()
        .function("set_replacement_model",            &modelFactory::setReplacementModel)
        .function("set_amino_replacement_model_file", &modelFactory::setCustomAAModelFile)
        .function("set_model_parameters",             &modelFactory::setModelParameters)
        .function("set_site_rate_model", optional_override([](modelFactory& self,
                const std::vector<MDOUBLE>& rates,
                const std::vector<MDOUBLE>& stationary_probs,
                const std::vector<std::vector<MDOUBLE>>& transition_matrix) {
            self.setSiteRateModel(rates, stationary_probs, transition_matrix);
        }))
        .function("set_site_rate_model_no_matrix", optional_override([](modelFactory& self,
                const std::vector<MDOUBLE>& rates,
                const std::vector<MDOUBLE>& stationary_probs) {
            self.setSiteRateModel(rates, stationary_probs, {});
        }))
        .function("reset",                &modelFactory::resetFactory)
        .function("is_model_valid",       &modelFactory::isModelValid)
        .function("build_replacement_model", &modelFactory::buildReplacementModel)
        .function("get_rate_category_sampler", optional_override([](modelFactory& self, MDOUBLE max_path_length) {
            return self.getRateCategorySampler(max_path_length);
        }), allow_raw_pointers())
        .function("get_rate_category_sampler_default", optional_override([](modelFactory& self) {
            return self.getRateCategorySampler(0);
        }), allow_raw_pointers());

    // -------------------------------------------------------------------------
    // SimulationContext
    //
    // protocol arg was optional (default nullptr) in pybind11.
    // Raw pointers require allow_raw_pointers() policy.
    // Two factory free-functions replace the overloaded constructor:
    //   createSimulationContext(tree, seed)
    //   createSimulationContextWithProtocol(tree, seed, protocol)
    //
    // JS caller is responsible for .delete() on the returned object.
    // -------------------------------------------------------------------------
    class_<SimulationContext<SelectedRNG>>("SimulationContext")
        .function("get_tree",            &SimulationContext<SelectedRNG>::getTree, allow_raw_pointers())
        .function("get_nodes_to_save",   &SimulationContext<SelectedRNG>::getNodesToSave)
        .function("set_save_leaves",     &SimulationContext<SelectedRNG>::setSaveLeaves)
        .function("set_save_root",       &SimulationContext<SelectedRNG>::setSaveRoot)
        .function("set_save_all",        &SimulationContext<SelectedRNG>::setSaveAll)
        .function("reseed",              &SimulationContext<SelectedRNG>::reseed)
        .function("set_branch_prob_cache", &SimulationContext<SelectedRNG>::setCacheBranchProbs)
        .function("get_indel_protocol",  &SimulationContext<SelectedRNG>::getProtocol, allow_raw_pointers())
        .function("set_protocol",        &SimulationContext<SelectedRNG>::setProtocol,  allow_raw_pointers())
        .function("set_category_sampler",&SimulationContext<SelectedRNG>::setCategorySampler, allow_raw_pointers())
        .function("get_category_sampler",&SimulationContext<SelectedRNG>::getCategorySampler, allow_raw_pointers());

    // Factory free-functions (replaces overloaded constructors)
    emscripten::function("createSimulationContext", optional_override([](tree* t, size_t seed) {
        return new SimulationContext<SelectedRNG>(t, seed, nullptr);
    }), allow_raw_pointers());

    emscripten::function("createSimulationContextWithProtocol", optional_override([](tree* t, size_t seed, SimulationProtocol* p) {
        return new SimulationContext<SelectedRNG>(t, seed, p);
    }), allow_raw_pointers());

    // -------------------------------------------------------------------------
    // IndelEvent / IndelEventType
    // -------------------------------------------------------------------------
    class_<Event>("IndelEvent")
        .property("type",     &Event::type)
        .property("position", &Event::position)
        .property("length",   &Event::length)
        .function("toString", optional_override([](const Event& e) {
            std::string type_name;
            switch (e.type) {
                case event::INSERTION: type_name = "INSERTION"; break;
                case event::DELETION:  type_name = "DELETION";  break;
                default:               type_name = "UNKNOWN";   break;
            }
            return "<IndelEvent type=" + type_name +
                   " position=" + std::to_string(e.position) +
                   " length="   + std::to_string(e.length) + ">";
        }));

    // -------------------------------------------------------------------------
    // IndelSimulator
    // -------------------------------------------------------------------------
    class_<IndelSimulator<SelectedRNG>>("IndelSimulator")
        .constructor<SimulationContext<SelectedRNG>&, SimulationProtocol*>(allow_raw_pointers())
        .function("update_protocol",  &IndelSimulator<SelectedRNG>::updateSimulationProtocol, allow_raw_pointers())
        .function("generate_events",  &IndelSimulator<SelectedRNG>::generateSimulation);

    // -------------------------------------------------------------------------
    // MSA
    // -------------------------------------------------------------------------
    class_<MSA<SelectedRNG>>("Msa")
        .function("length",                    &MSA<SelectedRNG>::getMSAlength)
        .function("num_sequences",             &MSA<SelectedRNG>::getNumberOfSequences)
        .function("fill_substitutions", optional_override([](MSA<SelectedRNG>& self, const SparseSequenceContainer& container) {
            self.fillSubstitutions(std::make_shared<const SparseSequenceContainer>(container));
        }))
        .function("print_msa",                 &MSA<SelectedRNG>::printFullMsa)
        .function("write_msa", optional_override([](MSA<SelectedRNG>& self, const std::string& path) {
            self.writeFullMsa(path.c_str());
        }))
        .function("get_msa_row_string",        &MSA<SelectedRNG>::generateMsaRowString)
        .function("get_sparse_msa",            &MSA<SelectedRNG>::getSparseMSA)
        .function("get_per_site_rate_categories", optional_override([](MSA<SelectedRNG>& self) {
            return *self.getPerSiteRateCategories();
        }))
        .function("get_root_positions_in_msa", &MSA<SelectedRNG>::getRootPositionsInMsa);

    // MSA factory free-functions (two constructors in pybind11)
    emscripten::function("createMsaFromEvents", optional_override([](EventMap& events, SimulationContext<SelectedRNG>& ctx) {
        return new MSA<SelectedRNG>(events, ctx);
    }), allow_raw_pointers());

    emscripten::function("createMsaFromLength", optional_override([](size_t length, SimulationContext<SelectedRNG>& ctx) {
        return new MSA<SelectedRNG>(length, ctx);
    }), allow_raw_pointers());

    // -------------------------------------------------------------------------
    // SubstitutionSimulator — Amino (20) and Nucleotide (4)
    // -------------------------------------------------------------------------
    using AminoSim = SubstitutionSimulator<SelectedRNG, 20>;
    class_<AminoSim>("AminoSubstitutionSimulator")
        .function("simulate_substitutions", optional_override([](AminoSim& self, size_t length) {
            return *self.simulateSubstitutions(length);
        }))
        .function("simulate_substitutions_with_root", optional_override([](AminoSim& self, size_t length, const std::string& rootSeq, const std::vector<size_t>& rootPositions) {
            return *self.simulateSubstitutions(length, rootSeq, rootPositions);
        }))
        .function("simulate_and_write_substitutions", &AminoSim::simulateAndWriteSubstitutions)
        .function("set_save_rates",                   &AminoSim::setSaveRates)
        .function("clear_rates_vec",                  &AminoSim::clearRatesVec)
        .function("set_aligned_sequence_map",         &AminoSim::setAlignedSequenceMap)
        .function("get_site_rates",                   &AminoSim::getSiteRates)
        .function("set_per_site_rate_categories", optional_override([](AminoSim& self, const std::vector<uint8_t>& cats) {
            self.setPerSiteRateCategories(std::make_shared<std::vector<uint8_t>>(cats));
        }))
        .function("get_per_site_rate_categories", optional_override([](AminoSim& self) {
            return *self.getPerSiteRateCategories();
        }));

    using NucleotideSim = SubstitutionSimulator<SelectedRNG, 4>;
    class_<NucleotideSim>("NucleotideSubstitutionSimulator")
        .function("simulate_substitutions", optional_override([](NucleotideSim& self, size_t length) {
            return *self.simulateSubstitutions(length);
        }))
        .function("simulate_substitutions_with_root", optional_override([](NucleotideSim& self, size_t length, const std::string& rootSeq, const std::vector<size_t>& rootPositions) {
            return *self.simulateSubstitutions(length, rootSeq, rootPositions);
        }))
        .function("simulate_and_write_substitutions", &NucleotideSim::simulateAndWriteSubstitutions)
        .function("set_save_rates",                   &NucleotideSim::setSaveRates)
        .function("clear_rates_vec",                  &NucleotideSim::clearRatesVec)
        .function("set_aligned_sequence_map",         &NucleotideSim::setAlignedSequenceMap)
        .function("get_site_rates",                   &NucleotideSim::getSiteRates)
        .function("set_per_site_rate_categories", optional_override([](NucleotideSim& self, const std::vector<uint8_t>& cats) {
            self.setPerSiteRateCategories(std::make_shared<std::vector<uint8_t>>(cats));
        }))
        .function("get_per_site_rate_categories", optional_override([](NucleotideSim& self) {
            return *self.getPerSiteRateCategories();
        }));

    emscripten::function("createAminoSim", optional_override([](modelFactory& factory, SimulationContext<SelectedRNG>& ctx) {
        return new AminoSim(factory, ctx);
    }), allow_raw_pointers());

    emscripten::function("createNucleotideSim", optional_override([](modelFactory& factory, SimulationContext<SelectedRNG>& ctx) {
        return new NucleotideSim(factory, ctx);
    }), allow_raw_pointers());
}