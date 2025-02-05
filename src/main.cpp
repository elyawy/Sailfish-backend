#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>

#include "./Simulator.h"

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

    py::class_<Block>(m, "Block")
        .def(py::init<size_t, size_t>());

    py::class_<BlockTree>(m, "BlockTree")
        .def(py::init<int>())
        .def("print_tree", &BlockTree::printTree)
        .def("block_list", &BlockTree::getBlockList);

    py::enum_<event>(m, "event")
        .value("Insertion", event::INSERTION)
        .value("Deletion", event::DELETION)
        .export_values();

    py::class_<DiscreteDistribution>(m, "DiscreteDistribution")
        .def(py::init<std::vector<double>>())
        .def("draw_sample", &DiscreteDistribution::drawSample, "Draw a random sample according to the given distribution")
        .def_static("set_seed", &DiscreteDistribution::setSeed, "Set seed for the random number generator")
        .def("get_table", &DiscreteDistribution::getTable, "Get Vose's alias table (useful for debugging)");

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
        .def(py::init<tree*>())
        .def("set_seed", &SimulationProtocol::setSeed)
        .def("get_seed", &SimulationProtocol::getSeed)
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

    py::enum_<alphabetCode>(m, "alphabetCode")
        .value("NULLCODE", alphabetCode::NULLCODE)
        .value("NUCLEOTIDE", alphabetCode::NUCLEOTIDE)
        .value("AMINOACID", alphabetCode::AMINOACID)
        .export_values();


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



    py::class_<modelFactory>(m, "modelFactory")
        .def(py::init<tree*>())
        .def("set_alphabet", &modelFactory::setAlphabet)
        .def("set_replacement_model" , &modelFactory::setReplacementModel)
        .def("set_amino_replacement_model_file" , &modelFactory::setCustomAAModelFile)
        .def("set_model_parameters" , &modelFactory::setModelParameters)
        .def("set_gamma_parameters" , &modelFactory::setGammaParameters)
        .def("set_invariant_sites_proportion", &modelFactory::setInvariantSitesProportion)
        .def("reset", &modelFactory::resetFactory);


    py::class_<Simulator>(m, "Simulator")
        .def(py::init<SimulationProtocol*>())
        .def("reset_sim", &Simulator::resetSimulator)
        .def("gen_indels", &Simulator::generateSimulation)
        .def("run_sim", &Simulator::runSimulator)
        .def("init_substitution_sim", &Simulator::initSubstitionSim)
        .def("gen_substitutions", &Simulator::simulateSubstitutions)
        .def("gen_substitutions_to_dir", &Simulator::simulateAndWriteSubstitutions)
        .def("save_site_rates", &Simulator::setSaveRates)
        .def("get_site_rates", &Simulator::getSiteRates)
        .def("save_all_nodes_sequences", &Simulator::setSaveAllNodes)
        .def("save_root_sequence", &Simulator::setSaveRoot)
        .def("get_saved_nodes_mask", &Simulator::getNodesSaveList);


    py::class_<MSA>(m, "Msa")
        .def(py::init<size_t, size_t, const std::vector<bool>& >())
        .def(py::init<BlockMap, tree::TreeNode*, const std::vector<bool>& >())
        .def("generate_msas", &MSA::generateMSAs)
        .def("length", &MSA::getMSAlength)
        .def("num_sequences", &MSA::getNumberOfSequences)
        .def("fill_substitutions", &MSA::fillSubstitutions)
        .def("set_substitutions_folder", &MSA::setSubstitutionsFolder)
        .def("print_msa", &MSA::printFullMsa)
        .def("print_indels", &MSA::printIndels)
        .def("write_msa", &MSA::writeFullMsa)
        .def("write_msa_from_dir", &MSA::writeMsaFromDir)
        .def("get_msa_string", &MSA::generateMsaString)
        .def("get_msa", &MSA::getMSAVec);

}
