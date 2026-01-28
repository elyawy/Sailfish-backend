#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "../../../libs/tomlplusplus/toml.hpp"
#include "../../../libs/pcg/pcg_random.hpp"
#include "../../../src/Simulator.h"
#include "../../../src/SimulationProtocol.h"
#include "../../../src/MSA.h"
#include "../../../src/modelFactory.h"
#include "../../../libs/Phylolib/includes/DiscreteDistribution.h"
#include "../../../libs/Phylolib/includes/tree.h"
#include "../../../libs/Phylolib/includes/gammaDistribution.h"

// Helper function to create a distribution from config
DiscreteDistribution* createDistribution(const toml::table& dist_config) {
    std::string type = dist_config["type"].value_or("zipf");
    double param = dist_config["parameter"].value_or(1.7);
    int truncation = dist_config["truncation"].value_or(50);
    
    std::vector<double> probs;
    
    if (type == "zipf") {
        // Zipf distribution: p(k) ∝ 1/k^param
        for (int i = 1; i <= truncation; ++i) {
            probs.push_back(1.0 / std::pow(i, param));
        }
    } else if (type == "geometric") {
        // Geometric distribution: p(k) = (1-p)^(k-1) * p
        for (int i = 1; i <= truncation; ++i) {
            probs.push_back(std::pow(1.0 - param, i - 1) * param);
        }
    } else if (type == "poisson") {
        // Poisson distribution: p(k) = (lambda^k * e^(-lambda)) / k!
        double lambda = param;
        double factorial = 1.0;
        for (int i = 1; i <= truncation; ++i) {
            factorial *= i;
            probs.push_back(std::pow(lambda, i) * std::exp(-lambda) / factorial);
        }
    } else {
        std::cerr << "Unknown distribution type: " << type << std::endl;
        exit(1);
    }
    
    // Normalize
    double sum = 0.0;
    for (auto p : probs) sum += p;
    for (auto& p : probs) p /= sum;
    
    return new DiscreteDistribution(probs);
}

// Helper to convert alphabet string to enum
alphabetCode getAlphabetCode(const std::string& alphabet) {
    if (alphabet == "protein") return alphabetCode::AMINOACID;
    if (alphabet == "nucleotide") return alphabetCode::NUCLEOTIDE;
    if (alphabet == "none") return alphabetCode::NULLCODE;
    
    std::cerr << "Unknown alphabet: " << alphabet << std::endl;
    exit(1);
}

// Helper to convert model string to enum
modelCode getModelCode(const std::string& model) {
    if (model == "NUCJC") return modelCode::NUCJC;
    if (model == "AAJC") return modelCode::AAJC;
    if (model == "GTR") return modelCode::GTR;
    if (model == "HKY") return modelCode::HKY;
    if (model == "TAMURA92") return modelCode::TAMURA92;
    if (model == "WAG") return modelCode::WAG;
    if (model == "LG") return modelCode::LG;
    if (model == "JONES") return modelCode::JONES;
    if (model == "DAYHOFF") return modelCode::DAYHOFF;
    if (model == "MTREV24") return modelCode::MTREV24;
    if (model == "CPREV45") return modelCode::CPREV45;
    if (model == "HIVB") return modelCode::HIVB;
    if (model == "HIVW") return modelCode::HIVW;
    
    std::cerr << "Unknown model: " << model << std::endl;
    exit(1);
}


int main(int argc, char* argv[]) {
    // Parse command line argument for config file
    std::string config_file = "test_config.toml";
    if (argc > 1) {
        config_file = argv[1];
    }
    
    std::cout << "Loading configuration from: " << config_file << std::endl;
    
    // Parse TOML config
    toml::table config;
    try {
        config = toml::parse_file(config_file);
    } catch (const toml::parse_error& err) {
        std::cerr << "Error parsing config file: " << err << std::endl;
        return 1;
    }
    
    // ==================== LOAD BASIC SETTINGS ====================
    std::string tree_file = config["basic"]["tree_file"].value_or("tree.nwk");
    int seed = config["basic"]["seed"].value_or(42);
    int root_seq_size = config["basic"]["root_sequence_size"].value_or(100);
    int min_seq_size = config["basic"]["min_sequence_size"].value_or(10);
    
    std::cout << "Tree file: " << tree_file << std::endl;
    std::cout << "Seed: " << seed << std::endl;
    std::cout << "Root sequence size: " << root_seq_size << std::endl;
    
    // ==================== LOAD TREE ====================
    tree phylo_tree(tree_file);
    size_t num_branches = phylo_tree.getNodesNum() - 1;
    std::cout << "Loaded tree with " << phylo_tree.getNodesNum() << " nodes" << std::endl;
    
    // ==================== SETUP PROTOCOL ====================
    SimulationProtocol protocol(&phylo_tree);
    protocol.setSeed(seed);
    protocol.setSequenceSize(root_seq_size);
    protocol.setMinSequenceSize(min_seq_size);
    
    // ==================== LOAD INDEL SETTINGS ====================
    bool no_indels = config["indels"]["no_indels"].value_or(false);
    double insertion_rate = config["indels"]["insertion_rate"].value_or(0.0);
    double deletion_rate = config["indels"]["deletion_rate"].value_or(0.0);
    
    std::cout << "Insertion rate: " << insertion_rate << std::endl;
    std::cout << "Deletion rate: " << deletion_rate << std::endl;
    
    // Set uniform rates for all branches
    std::vector<double> insertion_rates(num_branches, insertion_rate);
    std::vector<double> deletion_rates(num_branches, deletion_rate);
    protocol.setInsertionRates(insertion_rates);
    protocol.setDeletionRates(deletion_rates);
    
    // Create distributions
    auto insertion_dist_config = config["indels"]["insertion_distribution"];
    auto deletion_dist_config = config["indels"]["deletion_distribution"];
    
    DiscreteDistribution* ins_dist = createDistribution(*insertion_dist_config.as_table());
    DiscreteDistribution* del_dist = createDistribution(*deletion_dist_config.as_table());
    
    std::vector<DiscreteDistribution*> insertion_dists(num_branches, ins_dist);
    std::vector<DiscreteDistribution*> deletion_dists(num_branches, del_dist);
    protocol.setInsertionLengthDistributions(insertion_dists);
    protocol.setDeletionLengthDistributions(deletion_dists);
    
    // ==================== CREATE SIMULATOR ====================
    std::string alphabet_str = config["substitution"]["alphabet"].value_or("protein");
    alphabetCode alphabet = getAlphabetCode(alphabet_str);
    
    // Select appropriate simulator based on alphabet
    if (alphabet == alphabetCode::AMINOACID) {
        std::cout << "Creating protein simulator" << std::endl;
        
        Simulator<pcg64_fast, 20> sim(&protocol);
        
        // ==================== LOAD SUBSTITUTION MODEL ====================
        std::string model_str = config["substitution"]["model"].value_or("WAG");
        double gamma_alpha = config["substitution"]["gamma_alpha"].value_or(1.0);
        int gamma_categories = config["substitution"]["gamma_categories"].value_or(1);
        
        std::cout << "Substitution model: " << model_str << std::endl;
        std::cout << "Gamma alpha: " << gamma_alpha << ", categories: " << gamma_categories << std::endl;
        
        modelFactory mFac(&phylo_tree);
        mFac.setAlphabet(alphabet);
        mFac.setReplacementModel(getModelCode(model_str));
        
        // Create gamma distribution for site rates
        gammaDistribution gammaDist(gamma_alpha, gamma_categories);
        std::vector<double> rates;
        std::vector<double> probs;
        for (size_t i = 0; i < gammaDist.categories(); ++i) {
            rates.push_back(gammaDist.rates(i));
            probs.push_back(gammaDist.ratesProb(i));
        }
        
        // Set site rate model (empty transition matrix = independent rates)
        mFac.setSiteRateModel(rates, probs);
        
        if (!mFac.isModelValid()) {
            std::cerr << "Invalid model configuration!" << std::endl;
            return 1;
        }
        
        sim.initSubstitionSim(mFac);
        // set threshold for switching to Gillespie simulator
        double gillespieThreshold = config["substitution"]["gillespie_threshold"].value_or(1e-5);
        sim.setGillespieThreshold(gillespieThreshold);
        
        std::cout << "Set Gillespie threshold to: " << gillespieThreshold << std::endl;

        // ==================== OUTPUT SETTINGS ====================
        bool save_root = config["output"]["save_root_sequence"].value_or(false);
        bool save_all = config["output"]["save_all_nodes"].value_or(false);
        bool save_rates = config["output"]["save_site_rates"].value_or(false);
        
        if (save_root) {
            std::cout << "Saving root sequence" << std::endl;
            sim.setSaveRoot();
        }
        if (save_all) {
            std::cout << "Saving all node sequences" << std::endl;
            sim.setSaveAllNodes();
        }
        if (save_rates) {
            std::cout << "Saving site rates" << std::endl;
            sim.setSaveRates(true);
        }
        
        auto saveList = sim.getNodesSaveList();
        
        // ==================== RUN SIMULATION ====================
        std::cout << "\nRunning simulation..." << std::endl;
        BlockMap blockmap;
        MSA msa(0,0, std::vector<bool>()); // Dummy init
        int msa_length = config["basic"]["root_sequence_size"].value_or(100);
        if (no_indels) {
            std::cout << "Indels are disabled for this simulation." << std::endl;
        } else {
            blockmap = sim.generateSimulation();
            std::cout << "Generated indels" << std::endl;
            msa = MSA(blockmap, phylo_tree.getRoot(), saveList);
            msa_length = msa.getMSAlength();
        }
        
        std::cout << "MSA length: " << msa_length << std::endl;
        
        
        if (!no_indels) {
            auto seq_container = sim.simulateSubstitutions(msa_length);
            std::cout << "Generated substitutions" << std::endl;
            msa.fillSubstitutions(seq_container);
            sim.setAlignedSequenceMap(msa);
        } else {
            std::string output_file = config["output"]["output_file"].value_or("output.fasta");

            sim.simulateAndWriteSubstitutions(msa_length, output_file);
            std::cout << "Wrote substitutions to file: " << output_file << std::endl;
            return 0;
        }
        
        if (save_rates) {
            auto rates = sim.getSiteRates();
            std::cout << "\nSite rates (" << rates.size() << " sites):" << std::endl;
            for (size_t i = 0; i < std::min(rates.size(), size_t(20)); ++i) {
                std::cout << rates[i] << " ";
            }
            if (rates.size() > 20) std::cout << "...";
            std::cout << std::endl;
        }
        
        // ==================== OUTPUT ====================
        std::string output_mode = config["output"]["output_mode"].value_or("console");
        
        if (output_mode == "file") {
            std::string output_file = config["output"]["output_file"].value_or("output.fasta");
            std::cout << "\nWriting MSA to: " << output_file << std::endl;
            msa.writeFullMsa(output_file.c_str());
        } else {
            bool print_msa = config["output"]["print_msa"].value_or(true);
            bool print_indels = config["output"]["print_indels"].value_or(false);
            
            if (print_msa) {
                std::cout << "\n=== MSA ===" << std::endl;
                msa.printFullMsa();
            }
            if (print_indels) {
                std::cout << "\n=== Indels ===" << std::endl;
                msa.printIndels();
            }
        }
        
    } else if (alphabet == alphabetCode::NUCLEOTIDE) {
        std::cout << "Creating nucleotide simulator" << std::endl;
        
        Simulator<pcg64_fast, 4> sim(&protocol);
        
        // Similar setup for nucleotide models...
        std::string model_str = config["substitution"]["model"].value_or("NUCJC");
        double gamma_alpha = config["substitution"]["gamma_alpha"].value_or(1.0);
        int gamma_categories = config["substitution"]["gamma_categories"].value_or(1);
        
        modelFactory mFac(&phylo_tree);
        mFac.setAlphabet(alphabet);
        mFac.setReplacementModel(getModelCode(model_str));
        
        // Handle model parameters if provided (e.g., for HKY, GTR)
        if (config["substitution"]["model_parameters"].is_array()) {
            auto params_array = *config["substitution"]["model_parameters"].as_array();
            std::vector<double> params;
            for (auto& p : params_array) {
                params.push_back(p.value_or(0.0));
            }
            mFac.setModelParameters(params);
        }
        
        // Create gamma distribution for site rates
        gammaDistribution gammaDist(gamma_alpha, gamma_categories);
        std::vector<double> rates;
        std::vector<double> probs;
        for (size_t i = 0; i < gammaDist.categories(); ++i) {
            rates.push_back(gammaDist.rates(i));
            probs.push_back(gammaDist.ratesProb(i));
        }
        
        mFac.setSiteRateModel(rates, probs);
        
        if (!mFac.isModelValid()) {
            std::cerr << "Invalid model configuration!" << std::endl;
            return 1;
        }
        
        sim.initSubstitionSim(mFac);
        // set threshold for switching to Gillespie simulator
        double gillespieThreshold = config["substitution"]["gillespie_threshold"].value_or(1e-5);
        sim.setGillespieThreshold(gillespieThreshold);
        
        std::cout << "Set Gillespie threshold to: " << gillespieThreshold << std::endl;

        // ==================== OUTPUT SETTINGS ====================
        bool save_root = config["output"]["save_root_sequence"].value_or(false);
        bool save_all = config["output"]["save_all_nodes"].value_or(false);
        bool save_rates = config["output"]["save_site_rates"].value_or(false);
        
        if (save_root) {
            std::cout << "Saving root sequence" << std::endl;
            sim.setSaveRoot();
        }
        if (save_all) {
            std::cout << "Saving all node sequences" << std::endl;
            sim.setSaveAllNodes();
        }
        if (save_rates) {
            std::cout << "Saving site rates" << std::endl;
            sim.setSaveRates(true);
        }
        
        auto saveList = sim.getNodesSaveList();
        
        // ==================== RUN SIMULATION ====================
        std::cout << "\nRunning simulation..." << std::endl;
        BlockMap blockmap;
        MSA msa(0,0, std::vector<bool>()); // Dummy init
        int msa_length = config["basic"]["root_sequence_size"].value_or(100);
        if (no_indels) {
            std::cout << "Indels are disabled for this simulation." << std::endl;
        } else {
            blockmap = sim.generateSimulation();
            std::cout << "Generated indels" << std::endl;
            msa = MSA(blockmap, phylo_tree.getRoot(), saveList);
            msa_length = msa.getMSAlength();
        }
        
        std::cout << "MSA length: " << msa_length << std::endl;
        
        
        if (!no_indels) {
            auto seq_container = sim.simulateSubstitutions(msa_length);
            std::cout << "Generated substitutions" << std::endl;
            msa.fillSubstitutions(seq_container);
            sim.setAlignedSequenceMap(msa);
        } else {
            std::string output_file = config["output"]["output_file"].value_or("output.fasta");

            sim.simulateAndWriteSubstitutions(msa_length, output_file);
            std::cout << "Wrote substitutions to file: " << output_file << std::endl;
            // output rates if requested
            if (save_rates) {
                auto rates = sim.getSiteRates();
                std::cout << "\nSite rates (" << rates.size() << " sites):" << std::endl;
                for (size_t i = 0; i < std::min(rates.size(), size_t(20)); ++i) {
                    std::cout << rates[i] << " ";
                }
                if (rates.size() > 20) std::cout << "...";
                std::cout << std::endl;
            }
            return 0;
        }
        
        // Output
        std::string output_mode = config["output"]["output_mode"].value_or("console");
        if (output_mode == "file") {
            std::string output_file = config["output"]["output_file"].value_or("output.fasta");
            msa.writeFullMsa(output_file.c_str());
        } else {
            msa.printFullMsa();
        }
        if (save_rates) {
            auto rates = sim.getSiteRates();
            std::cout << "\nSite rates (" << rates.size() << " sites):" << std::endl;
            for (size_t i = 0; i < std::min(rates.size(), size_t(20)); ++i) {
                std::cout << rates[i] << " ";
            }
            if (rates.size() > 20) std::cout << "...";
            std::cout << std::endl;
        }
        
    } else {
        // Indels only (no substitutions)
        std::cout << "Running indels-only simulation" << std::endl;
        
        Simulator<pcg64_fast, 20> sim(&protocol);
        
        bool save_root = config["output"]["save_root_sequence"].value_or(false);
        if (save_root) sim.setSaveRoot();
        
        auto saveList = sim.getNodesSaveList();
        
        auto blockmap = sim.generateSimulation();
        MSA msa(blockmap, phylo_tree.getRoot(), saveList);
        
        std::string output_mode = config["output"]["output_mode"].value_or("console");
        if (output_mode == "file") {
            std::string output_file = config["output"]["output_file"].value_or("output.fasta");
            std::cout << "Writing indels to: " << output_file << std::endl;
            // Note: For indels-only, we'd need a different output method
            msa.printIndels();
        } else {
            msa.printIndels();
        }
    }
    
    std::cout << "\nSimulation complete!" << std::endl;
    
    // Cleanup
    delete ins_dist;
    delete del_dist;
    
    return 0;
}