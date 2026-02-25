#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <memory>

#include "../libs/pcg/pcg_random.hpp"
#include "../libs/Phylolib/includes/gammaDistribution.h"
#include "./IndelSimulator.h"
#include "./SubstitutionSimulator.h"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(SparseSequenceContainer);


PYBIND11_MODULE(_Sailfish, m) {
    m.doc() = R"pbdoc(
        Sailfish simulator
        -----------------------

        .. currentmodule:: _Sailfish

        .. autosummary::
           :toctree: _generate

            DiscreteDistribution
            SimProtocol
            alphabetCode
            modelCode
            modelFactory
            Simulator
            Msa
            Tree
    )pbdoc";

    using SelectedRNG = pcg64_fast;

    py::class_<DiscreteDistribution>(m, "DiscreteDistribution")
        .def(py::init<std::vector<double>>());

    py::bind_vector<SparseSequenceContainer, std::shared_ptr<SparseSequenceContainer>>(m, "SparseSequenceContainer");


    py::class_<tree>(m, "Tree")
        .def(py::init<const std::string&, bool>(), "Create Phylogenetic tree object from newick formatted file")
        .def_property_readonly("num_nodes", &tree::getNodesNum)
        .def_property_readonly("root", &tree::getRoot);

    py::class_<tree::TreeNode>(m, "node")
        .def_property_readonly("sons", &tree::TreeNode::getSons)
        .def_property_readonly("num_leaves", &tree::TreeNode::getNumberLeaves)
        .def_property_readonly("name", &tree::TreeNode::name)
        .def("distance_to_father", &tree::TreeNode::dis2father);


    py::enum_<SiteRateModel>(m, "SiteRateModel")
        .value("SIMPLE", SiteRateModel::SIMPLE)
        .value("INDEL_AWARE", SiteRateModel::INDEL_AWARE)
        .export_values();

    py::class_<SimulationProtocol>(m, "SimProtocol")
        .def(py::init<size_t>())
        .def("set_sequence_size", &SimulationProtocol::setSequenceSize)
        .def("get_sequence_size", &SimulationProtocol::getSequenceSize)
        .def("set_insertion_rates", &SimulationProtocol::setInsertionRates)
        .def("get_insertion_rate", &SimulationProtocol::getInsertionRate)
        .def("set_deletion_rates", &SimulationProtocol::setDeletionRates)
        .def("get_deletion_rate", &SimulationProtocol::getDeletionRate)
        .def("set_insertion_length_distributions", &SimulationProtocol::setInsertionLengthDistributions)
        .def("get_insertion_length_distribution", &SimulationProtocol::getInsertionDistribution)
        .def("set_deletion_length_distributions", &SimulationProtocol::setDeletionLengthDistributions)
        .def("get_deletion_length_distribution", &SimulationProtocol::getDeletionDistribution)
        .def("set_minimum_sequence_size", &SimulationProtocol::setMinSequenceSize)
        .def("get_minimum_sequence_size", &SimulationProtocol::getMinSequenceSize)
        .def("set_site_rate_model", &SimulationProtocol::setSiteRateModel)
        .def("get_site_rate_model", &SimulationProtocol::getSiteRateModel)
        .def("set_max_insertion_length", &SimulationProtocol::setMaxInsertionLength)
        .def("get_max_insertion_length", &SimulationProtocol::getMaxInsertionLength);



    // py::class_<sequenceContainer, std::shared_ptr<sequenceContainer>>(m, "sequenceContainer")
    //     .def(py::init<>());


    py::enum_<modelCode>(m, "modelCode")
        .value("NUCJC", modelCode::NUCJC)
        .value("AAJC", modelCode::AAJC)
        .value("GTR", modelCode::GTR)
        .value("HKY", modelCode::HKY)
        .value("TAMURA92", modelCode::TAMURA92)
        // .value("WYANGMODEL", modelCode::WYANGMODEL)
        .value("CPREV45", modelCode::CPREV45)
        .value("DAYHOFF", modelCode::DAYHOFF)
        .value("JONES", modelCode::JONES)	// THIS IS JTT
        .value("MTREV24", modelCode::MTREV24)
        .value("WAG", modelCode::WAG)
        .value("HIVB", modelCode::HIVB)
        .value("HIVW", modelCode::HIVW)
        .value("LG", modelCode::LG)
        .value("EMPIRICODON", modelCode::EMPIRICODON)
        .value("EX_BURIED", modelCode::EX_BURIED) 
        .value("EX_EXPOSED", modelCode::EX_EXPOSED)
        .value("EHO_EXTENDED", modelCode::EHO_EXTENDED)
        .value("EHO_HELIX", modelCode::EHO_HELIX)
        .value("EHO_OTHER", modelCode::EHO_OTHER)
        .value("EX_EHO_BUR_EXT", modelCode::EX_EHO_BUR_EXT)
        .value("EX_EHO_BUR_HEL", modelCode::EX_EHO_BUR_HEL)
        .value("EX_EHO_BUR_OTH", modelCode::EX_EHO_BUR_OTH)
        .value("EX_EHO_EXP_EXT", modelCode::EX_EHO_EXP_EXT)
        .value("EX_EHO_EXP_HEL", modelCode::EX_EHO_EXP_HEL)
        .value("EX_EHO_EXP_OTH", modelCode::EX_EHO_EXP_OTH)
        .value("CUSTOM", modelCode::CUSTOM)
        .export_values();

    py::class_<gammaDistribution>(m, "GammaDistribution")
        .def(py::init<MDOUBLE, int>())
        .def("getAllRates", [](const gammaDistribution& g) {
            std::vector<MDOUBLE> result;
            for (size_t i = 0; i < g.categories(); ++i) {
                result.push_back(g.rates(i));
            }
            return result;
        })
        .def("getAllRatesProb", [](const gammaDistribution& g) {
            std::vector<MDOUBLE> result;
            for (size_t i = 0; i < g.categories(); ++i) {
                result.push_back(g.ratesProb(i));
            }
            return result;
        });

    py::class_<CategorySampler>(m, "CategorySampler");

    py::class_<modelFactory>(m, "modelFactory")
        .def(py::init<>())
        .def("set_replacement_model" , &modelFactory::setReplacementModel)
        .def("set_amino_replacement_model_file" , &modelFactory::setCustomAAModelFile)
        .def("set_model_parameters" , &modelFactory::setModelParameters)
        .def("set_site_rate_model", &modelFactory::setSiteRateModel,
            py::arg("rates"),
            py::arg("stationary_probs"),
            py::arg("transition_matrix") = std::vector<std::vector<MDOUBLE>>())
        .def("reset", &modelFactory::resetFactory)
        .def("is_model_valid", &modelFactory::isModelValid)
        .def("build_replacement_model", &modelFactory::buildReplacementModel)
        .def("get_rate_category_sampler", &modelFactory::getRateCategorySampler, py::arg("max_path_length") = 0);

    py::class_<SimulationContext<SelectedRNG>>(m, "SimulationContext")
        .def(py::init<tree*, size_t, SimulationProtocol*>(),
            py::arg("tree"), py::arg("seed"), 
            py::arg("protocol") = nullptr)
        .def("get_tree", &SimulationContext<SelectedRNG>::getTree)
        .def("get_nodes_to_save", &SimulationContext<SelectedRNG>::getNodesToSave)
        .def("set_save_leaves", &SimulationContext<SelectedRNG>::setSaveLeaves)
        .def("set_save_root", &SimulationContext<SelectedRNG>::setSaveRoot)
        .def("set_save_all", &SimulationContext<SelectedRNG>::setSaveAll)
        .def("reseed", &SimulationContext<SelectedRNG>::reseed)
        .def("get_indel_protocol", &SimulationContext<SelectedRNG>::getProtocol)
        .def("set_protocol", &SimulationContext<SelectedRNG>::setProtocol)
        .def("set_category_sampler", &SimulationContext<SelectedRNG>::setCategorySampler)
        .def("get_category_sampler", &SimulationContext<SelectedRNG>::getCategorySampler);

    //event enum bindings
    py::enum_<event>(m, "IndelEventType")
        .value("INSERTION", event::INSERTION)
        .value("DELETION", event::DELETION)
        .export_values();

    //Indel event struct bindings
    py::class_<Event>(m, "IndelEvent")
        .def_readonly("type", &Event::type)
        .def_readonly("position", &Event::position)
        .def_readonly("length", &Event::length)
        .def("__repr__", [](const Event& e) {
            std::string type_name;
            switch(e.type) {
                case event::INSERTION: type_name = "INSERTION"; break;
                case event::DELETION: type_name = "DELETION"; break;
                default: type_name = "UNKNOWN"; break;
            }
            return "<IndelEvent type=" + type_name + 
                " position=" + std::to_string(e.position) + 
                " length=" + std::to_string(e.length) + ">";
        });

    // bindings for IndelSimulator
    py::class_<IndelSimulator<SelectedRNG>>(m, "IndelSimulator")
        .def(py::init<SimulationContext<SelectedRNG>&, SimulationProtocol*>())
        .def("update_protocol", &IndelSimulator<SelectedRNG>::updateSimulationProtocol)
        .def("generate_events", &IndelSimulator<SelectedRNG>::generateSimulation);

    // bindings for SubstitutionSimulator (amino)
    py::class_<SubstitutionSimulator<SelectedRNG, 20>>(m, "AminoSubstitutionSimulator")
        .def(py::init<modelFactory&, SimulationContext<SelectedRNG>&>())
        .def("simulate_substitutions", [](SubstitutionSimulator<SelectedRNG, 20>& self, size_t length) {
            return *self.simulateSubstitutions(length);})
        .def("simulate_and_write_substitutions", &SubstitutionSimulator<SelectedRNG, 20>::simulateAndWriteSubstitutions)
        .def("init_substitution_sim", &SubstitutionSimulator<SelectedRNG, 20>::initSubstitionSim)
        .def("set_save_rates", &SubstitutionSimulator<SelectedRNG, 20>::setSaveRates)
        .def("clear_rates_vec", &SubstitutionSimulator<SelectedRNG, 20>::clearRatesVec)
        // .def("get_sequence_container", &SubstitutionSimulator<SelectedRNG, 20>::getSequenceContainer)
        .def("set_aligned_sequence_map", &SubstitutionSimulator<SelectedRNG, 20>::setAlignedSequenceMap)
        .def("get_site_rates", &SubstitutionSimulator<SelectedRNG, 20>::getSiteRates)
        .def("set_per_site_rate_categories", [](SubstitutionSimulator<SelectedRNG, 20>& self, std::vector<size_t> cats) {
            self.setPerSiteRateCategories(std::make_shared<const std::vector<size_t>>(std::move(cats)));
        })
        .def("get_per_site_rate_categories", [](SubstitutionSimulator<SelectedRNG, 20>& self) {
            auto cats = self.getPerSiteRateCategories();
            return cats ? std::vector<size_t>(*cats) : std::vector<size_t>{};
        });

    // bindings for SubstitutionSimulator (nucleotide)
    py::class_<SubstitutionSimulator<SelectedRNG, 4>>(m, "NucleotideSubstitutionSimulator")
        .def(py::init<modelFactory&, SimulationContext<SelectedRNG>&>())
        .def("simulate_substitutions", [](SubstitutionSimulator<SelectedRNG, 4>& self, size_t length) {
            return *self.simulateSubstitutions(length);})
        .def("simulate_and_write_substitutions", &SubstitutionSimulator<SelectedRNG, 4>::simulateAndWriteSubstitutions)
        .def("init_substitution_sim", &SubstitutionSimulator<SelectedRNG, 4>::initSubstitionSim)
        .def("set_save_rates", &SubstitutionSimulator<SelectedRNG, 4>::setSaveRates)
        .def("clear_rates_vec", &SubstitutionSimulator<SelectedRNG, 4>::clearRatesVec)
        // .def("get_sequence_container", &SubstitutionSimulator<SelectedRNG, 4>::getSequenceContainer)
        .def("set_aligned_sequence_map", &SubstitutionSimulator<SelectedRNG, 4>::setAlignedSequenceMap)
        .def("get_site_rates", &SubstitutionSimulator<SelectedRNG, 4>::getSiteRates)
        .def("set_per_site_rate_categories", [](SubstitutionSimulator<SelectedRNG, 4>& self, std::vector<size_t> cats) {
            self.setPerSiteRateCategories(std::make_shared<const std::vector<size_t>>(std::move(cats)));
        })
        .def("get_per_site_rate_categories", [](SubstitutionSimulator<SelectedRNG, 4>& self) {
            auto cats = self.getPerSiteRateCategories();
            return cats ? std::vector<size_t>(*cats) : std::vector<size_t>{};
        });


    py::class_<MSA<SelectedRNG>>(m, "Msa")
        .def(py::init<EventMap&, SimulationContext<SelectedRNG>&>())
        .def(py::init<size_t, SimulationContext<SelectedRNG>&>())
        .def("length", &MSA<SelectedRNG>::getMSAlength)
        .def("num_sequences", &MSA<SelectedRNG>::getNumberOfSequences)
        .def("fill_substitutions", &MSA<SelectedRNG>::fillSubstitutions)
        .def("print_msa", &MSA<SelectedRNG>::printFullMsa)
        .def("write_msa", &MSA<SelectedRNG>::writeFullMsa)
        .def("get_msa_row_string", &MSA<SelectedRNG>::generateMsaRowString)
        .def("get_sparse_msa", &MSA<SelectedRNG>::getSparseMSA)
        .def("get_per_site_rate_categories", [](MSA<SelectedRNG>& self) {
            auto cats = self.getPerSiteRateCategories();
            return cats ? *cats : std::vector<size_t>{};
        });

}
