#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>

#include "../libs/pcg/pcg_random.hpp"
#include "../libs/Phylolib/includes/gammaDistribution.h"
#include "./IndelSimulator.h"
#include "./SubstitutionSimulator.h"

namespace py = pybind11;


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

    py::class_<tree>(m, "Tree")
        .def(py::init<const std::string&, bool>(), "Create Phylogenetic tree object from newick formatted file")
        .def_property_readonly("num_nodes", &tree::getNodesNum)
        .def_property_readonly("root", &tree::getRoot);

    py::class_<tree::TreeNode>(m, "node")
        .def_property_readonly("sons", &tree::TreeNode::getSons)
        .def_property_readonly("num_leaves", &tree::TreeNode::getNumberLeaves)
        .def_property_readonly("name", &tree::TreeNode::name)
        .def("distance_to_father", &tree::TreeNode::dis2father);

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
        .def("get_minimum_sequence_size", &SimulationProtocol::getMinSequenceSize);



    py::class_<sequenceContainer, std::shared_ptr<sequenceContainer>>(m, "sequenceContainer")
        .def(py::init<>());


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

    py::class_<modelFactory>(m, "modelFactory")
        .def(py::init<>())
        .def("set_replacement_model" , &modelFactory::setReplacementModel)
        .def("set_amino_replacement_model_file" , &modelFactory::setCustomAAModelFile)
        .def("set_model_parameters" , &modelFactory::setModelParameters)
        .def("setSiteRateModel", &modelFactory::setSiteRateModel,
            py::arg("rates"),
            py::arg("stationary_probs"),
            py::arg("transition_matrix") = std::vector<std::vector<MDOUBLE>>())
        .def("reset", &modelFactory::resetFactory);


    py::class_<SimulationContext<SelectedRNG>>(m, "SimulationContext")
        .def(py::init<tree*, size_t>())
        .def("get_tree", &SimulationContext<SelectedRNG>::getTree)
        .def("get_nodes_to_save", &SimulationContext<SelectedRNG>::getNodesToSave);

    // bindings for IndelSimulator
    py::class_<IndelSimulator<SelectedRNG>>(m, "IndelSimulator")
        .def(py::init<SimulationContext<SelectedRNG>&, SimulationProtocol*>())
        .def("update_protocol", &IndelSimulator<SelectedRNG>::updateSimulationProtocol)
        .def("generate_simulation", &IndelSimulator<SelectedRNG>::generateSimulation);

    // bindings for SubstitutionSimulator (amino)
    py::class_<SubstitutionSimulator<SelectedRNG, 20>>(m, "AminoSubstitutionSimulator")
        .def(py::init<modelFactory&, SimulationContext<SelectedRNG>&>())
        .def("init_substitution_sim", &SubstitutionSimulator<SelectedRNG, 20>::initSubstitionSim)
        .def("set_save_rates", &SubstitutionSimulator<SelectedRNG, 20>::setSaveRates)
        .def("clear_rates_vec", &SubstitutionSimulator<SelectedRNG, 20>::clearRatesVec)
        .def("get_sequence_container", &SubstitutionSimulator<SelectedRNG, 20>::getSequenceContainer)
        .def("get_site_rates", &SubstitutionSimulator<SelectedRNG, 20>::getSiteRates);

    // bindings for SubstitutionSimulator (nucleotide)
    py::class_<SubstitutionSimulator<SelectedRNG, 4>>(m, "NucleotideSubstitutionSimulator")
        .def(py::init<modelFactory&, SimulationContext<SelectedRNG>&>())
        .def("init_substitution_sim", &SubstitutionSimulator<SelectedRNG, 4>::initSubstitionSim)
        .def("set_save_rates", &SubstitutionSimulator<SelectedRNG, 4>::setSaveRates)
        .def("clear_rates_vec", &SubstitutionSimulator<SelectedRNG, 4>::clearRatesVec)
        .def("get_sequence_container", &SubstitutionSimulator<SelectedRNG, 4>::getSequenceContainer)
        .def("get_site_rates", &SubstitutionSimulator<SelectedRNG, 4>::getSiteRates);


    py::class_<MSA>(m, "Msa")
        .def(py::init<size_t, size_t, const std::vector<bool>& >())
        .def(py::init<EventMap&, size_t, tree::TreeNode*, const std::vector<bool>& >())
        .def("length", &MSA::getMSAlength)
        .def("num_sequences", &MSA::getNumberOfSequences)
        .def("fill_substitutions", &MSA::fillSubstitutions)
        .def("print_msa", &MSA::printFullMsa)
        .def("write_msa", &MSA::writeFullMsa)
        .def("get_msa_string", &MSA::generateMsaString)
        .def("get_msa", &MSA::getMSAVec);

}
