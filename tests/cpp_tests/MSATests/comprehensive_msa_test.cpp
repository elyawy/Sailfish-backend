#include <random>
#include <sstream>
#include <fstream>
#include "../../../src/Simulator.h"
#include "../../../src/MSA.h"
#include "../../../src/MsaFixed.h"

// Helper to generate random Newick tree
std::string generateRandomTree(std::mt19937& rng, int numLeaves) {
    if (numLeaves == 1) {
        static int leafId = 0;
        return "L" + std::to_string(leafId++);
    }
    
    std::uniform_int_distribution<> splitDist(1, numLeaves - 1);
    std::uniform_real_distribution<> branchDist(0.01, 0.5);
    
    int leftSize = splitDist(rng);
    int rightSize = numLeaves - leftSize;
    
    std::stringstream ss;
    ss << "(" << generateRandomTree(rng, leftSize) << ":" << branchDist(rng)
       << "," << generateRandomTree(rng, rightSize) << ":" << branchDist(rng) << ")";
    
    return ss.str();
}

// Reset static counter between trials
void resetLeafCounter() {
    static int dummy = 0;
    dummy = 0;
}

// Generate a standalone test program for the failing case
void generateFailingTestProgram(const std::string& newick, int rootLength,
                                double insertionRate, double deletionRate,
                                int simSeed, int trialNumber) {
    std::ofstream out("failing_test_trial_" + std::to_string(trialNumber) + ".cpp");
    
    out << "#include \"../../../src/Simulator.h\"\n";
    out << "#include \"../../../src/MSA.h\"\n";
    out << "#include \"../../../src/MsaFixed.h\"\n\n";
    out << "int main() {\n";
    out << "    // Failed trial " << trialNumber << "\n";
    out << "    std::string newick = \"" << newick << "\";\n";
    out << "    int rootLength = " << rootLength << ";\n";
    out << "    double insertionRate = " << insertionRate << ";\n";
    out << "    double deletionRate = " << deletionRate << ";\n";
    out << "    int simSeed = " << simSeed << ";\n\n";
    
    out << "    tree tree_(newick, false);\n\n";
    
    out << "    vector<DiscreteDistribution*> insertionDists(tree_.getNodesNum() - 1);\n";
    out << "    vector<DiscreteDistribution*> deletionDists(tree_.getNodesNum() - 1);\n\n";
    
    out << "    DiscreteDistribution indelDist({0.500597013622622,0.15407680416466002,0.07733607404815325,0.04742269916833997,0.03245183409648796,0.023802968878603163,0.018315618393686197,0.014596047786709592,0.011947471092366838,0.009988223562671806,0.008494164935846923,0.007326223038978798,0.00639417555007358,0.005637292815586038,0.005013408742738382,0.004492460672380736,0.004052519265065782,0.0036772655722418293,0.0033543403993036617,0.0030742364096058624,0.00282953749580224,0.002614385926699555,0.0024241022418638207,0.0022549096413393907,0.0021037311601300934,0.001968038376604936,0.0018457371577210268,0.001735080389074546,0.001634600615111441,0.00154305754135923,0.001459396750517745,0.0013827169647433223,0.0013122438819642994,0.0012473091132895959,0.0011873331108893934,0.0011318112414127814,0.0010803023567287186,0.0010324193607649647,0.0009878213819912308,0.000946207245246948,0.0009073100010352715,0.000870892320076723,0.0008367425994772478,0.0008046716569959321,0.0007745099135856162,0.0007461049831118826,0.0007193196030519133,0.000694029851883145,0.0006701236084401137,0.000647499216243446,0.0006260643220715613,0.0006057348631513914,0.0005864341815220788,0.0005680922475578417,0.0005506449774670669,0.0005340336319283628,0.0005182042849716617,0.0005031073538361537,0.0004886971818952107,0.00047493166787853947,0.0004617719355816668,0.00044918203906341153,0.00043712869901847034,0.00042558106659524523,0.0004145105114255328,0.0004038904310565982,0.0003936960793390335,0.00038390441163517177,0.00037449394498066375,0.0003654446315627292,0.0003567377440781168,0.00034835577170658,0.00034028232558560196,0.0003325020528024772,0.00032500055803344507,0.000317764332058729,0.0003107806864690648,0.0003040376939552896,0.00029752413363926596,0.00029122944096306576,0.0002851436617049981,0.00027925740973663626,0.00027356182817527154,0.0002680485536218626,0.00026270968320613465,0.0002575377441885291,0.00025252566589363025,0.00024766675377188466,0.00024295466540621684,0.00023838338829779629,0.00023394721928100186,0.00022964074543174937,0.0002254588263460113,0.00022139657767671146,0.00021744935582738057,0.00021361274371013342,0.00020988253748378816,0.00020625473419539554,0.00020272552025516334,0.0001992912606808326});\n\n";
    
    out << "    fill(insertionDists.begin(), insertionDists.end(), &indelDist);\n";
    out << "    fill(deletionDists.begin(), deletionDists.end(), &indelDist);\n\n";
    
    out << "    vector<double> insertionRates(tree_.getNodesNum() - 1);\n";
    out << "    vector<double> deletionRates(tree_.getNodesNum() - 1);\n";
    out << "    fill(insertionRates.begin(), insertionRates.end(), insertionRate);\n";
    out << "    fill(deletionRates.begin(), deletionRates.end(), deletionRate);\n\n";
    
    out << "    SimulationProtocol protocol(&tree_);\n";
    out << "    protocol.setInsertionLengthDistributions(insertionDists);\n";
    out << "    protocol.setDeletionLengthDistributions(deletionDists);\n";
    out << "    protocol.setInsertionRates(insertionRates);\n";
    out << "    protocol.setDeletionRates(deletionRates);\n";
    out << "    protocol.setSequenceSize(rootLength);\n";
    out << "    protocol.setSeed(simSeed);\n\n";
    
    out << "    Simulator sim(&protocol);\n";
    out << "    auto saveList = sim.getNodesSaveList();\n";
    out << "    auto blockmap = sim.generateSimulation();\n\n";
    
    out << "    // Build with both implementations\n";
    out << "    MSA msaOriginal(blockmap, tree_.getRoot(), saveList);\n";
    out << "    MsaFixed msaFixed(blockmap, tree_.getRoot(), saveList);\n\n";
    
    out << "    std::string strOriginal = msaOriginal.generateMsaString();\n";
    out << "    std::string strFixed = msaFixed.generateMsaString();\n\n";
    
    out << "    std::cout << \"Original MSA length: \" << msaOriginal.getMSAlength() << \"\\n\";\n";
    out << "    std::cout << \"Fixed MSA length: \" << msaFixed.getMSAlength() << \"\\n\\n\";\n\n";
    
    out << "    if (strOriginal != strFixed) {\n";
    out << "        std::cout << \"✗ MISMATCH REPRODUCED!\\n\";\n";
    out << "        std::cout << \"\\n--- Original MSA ---\\n\" << strOriginal << \"\\n\";\n";
    out << "        std::cout << \"\\n--- Fixed MSA ---\\n\" << strFixed << \"\\n\";\n";
    out << "        return 1;\n";
    out << "    } else {\n";
    out << "        std::cout << \"✓ No mismatch found (may be environment-dependent)\\n\";\n";
    out << "        return 0;\n";
    out << "    }\n";
    out << "}\n";
    
    out.close();
    
    std::cout << "\n*** Generated standalone test: failing_test_trial_" << trialNumber << ".cpp ***\n";
}

int main() {
    std::random_device rd;
    std::mt19937 rng(rd());
    
    // Distribution for varying parameters
    std::uniform_int_distribution<> numLeavesDist(2, 10);
    std::uniform_int_distribution<> rootLengthDist(10, 100);
    std::uniform_real_distribution<> indelRateDist(0.0, 1.0);
    std::uniform_int_distribution<> seedDist(0, 999999);
    
    // Zipf distribution with 'a' parameter of 1.7
    DiscreteDistribution indelDist({0.500597013622622,0.15407680416466002,0.07733607404815325,0.04742269916833997,0.03245183409648796,0.023802968878603163,0.018315618393686197,0.014596047786709592,0.011947471092366838,0.009988223562671806,0.008494164935846923,0.007326223038978798,0.00639417555007358,0.005637292815586038,0.005013408742738382,0.004492460672380736,0.004052519265065782,0.0036772655722418293,0.0033543403993036617,0.0030742364096058624,0.00282953749580224,0.002614385926699555,0.0024241022418638207,0.0022549096413393907,0.0021037311601300934,0.001968038376604936,0.0018457371577210268,0.001735080389074546,0.001634600615111441,0.00154305754135923,0.001459396750517745,0.0013827169647433223,0.0013122438819642994,0.0012473091132895959,0.0011873331108893934,0.0011318112414127814,0.0010803023567287186,0.0010324193607649647,0.0009878213819912308,0.000946207245246948,0.0009073100010352715,0.000870892320076723,0.0008367425994772478,0.0008046716569959321,0.0007745099135856162,0.0007461049831118826,0.0007193196030519133,0.000694029851883145,0.0006701236084401137,0.000647499216243446,0.0006260643220715613,0.0006057348631513914,0.0005864341815220788,0.0005680922475578417,0.0005506449774670669,0.0005340336319283628,0.0005182042849716617,0.0005031073538361537,0.0004886971818952107,0.00047493166787853947,0.0004617719355816668,0.00044918203906341153,0.00043712869901847034,0.00042558106659524523,0.0004145105114255328,0.0004038904310565982,0.0003936960793390335,0.00038390441163517177,0.00037449394498066375,0.0003654446315627292,0.0003567377440781168,0.00034835577170658,0.00034028232558560196,0.0003325020528024772,0.00032500055803344507,0.000317764332058729,0.0003107806864690648,0.0003040376939552896,0.00029752413363926596,0.00029122944096306576,0.0002851436617049981,0.00027925740973663626,0.00027356182817527154,0.0002680485536218626,0.00026270968320613465,0.0002575377441885291,0.00025252566589363025,0.00024766675377188466,0.00024295466540621684,0.00023838338829779629,0.00023394721928100186,0.00022964074543174937,0.0002254588263460113,0.00022139657767671146,0.00021744935582738057,0.00021361274371013342,0.00020988253748378816,0.00020625473419539554,0.00020272552025516334,0.0001992912606808326});

    int numTrials = 1000;
    int numMismatches = 0;
    
    for (int trial = 0; trial < numTrials; trial++) {
        resetLeafCounter();
        
        // Generate random parameters
        int numLeaves = numLeavesDist(rng);
        int rootLength = rootLengthDist(rng);
        double insertionRate = indelRateDist(rng);
        double deletionRate = indelRateDist(rng);
        int simSeed = seedDist(rng);
        
        std::string newick = generateRandomTree(rng, numLeaves) + ";";
        
        std::cout << "\n=== Trial " << (trial + 1) << "/" << numTrials << " ===\n";
        std::cout << "Leaves: " << numLeaves << ", Root length: " << rootLength 
                  << ", Ins rate: " << insertionRate << ", Del rate: " << deletionRate
                  << ", Seed: " << simSeed << "\n";
        std::cout << "Tree: " << newick << "\n";
        
        try {
            tree tree_(newick, false);
            
            vector<DiscreteDistribution*> insertionDists(tree_.getNodesNum() - 1);
            vector<DiscreteDistribution*> deletionDists(tree_.getNodesNum() - 1);
            
            fill(insertionDists.begin(), insertionDists.end(), &indelDist);
            fill(deletionDists.begin(), deletionDists.end(), &indelDist);
            
            vector<double> insertionRates(tree_.getNodesNum() - 1);
            vector<double> deletionRates(tree_.getNodesNum() - 1);
            fill(insertionRates.begin(), insertionRates.end(), insertionRate);
            fill(deletionRates.begin(), deletionRates.end(), deletionRate);
            
            SimulationProtocol protocol(&tree_);
            protocol.setInsertionLengthDistributions(insertionDists);
            protocol.setDeletionLengthDistributions(deletionDists);
            protocol.setInsertionRates(insertionRates);
            protocol.setDeletionRates(deletionRates);
            protocol.setSequenceSize(rootLength);
            protocol.setSeed(simSeed);
            
            Simulator sim(&protocol);
            auto saveList = sim.getNodesSaveList();
            
            auto blockmap = sim.generateSimulation();
            // Build with both implementations
            
            MSA msaOriginal(blockmap, tree_.getRoot(), saveList);

            MsaFixed msaFixed(blockmap, tree_.getRoot(), saveList);

            std::string strOriginal = msaOriginal.generateMsaString();
            std::string strFixed = msaFixed.generateMsaString();
            
            // Compare
            if (strOriginal != strFixed) {
                numMismatches++;
                std::cout << "\n✗ MISMATCH FOUND!\n";
                std::cout << "Original MSA length: " << msaOriginal.getMSAlength() << "\n";
                std::cout << "Fixed MSA length: " << msaFixed.getMSAlength() << "\n";
                
                // Generate standalone test program
                generateFailingTestProgram(newick, rootLength, insertionRate, 
                                          deletionRate, simSeed, trial + 1);
                
                // Print blockmaps
                std::cout << "\n--- Block Maps ---\n";
                for (const auto& [nodeId, blockTuple] : blockmap) {
                    const BlockList& blocks = std::get<static_cast<int>(BLOCKLIST::BLOCKS)>(blockTuple);
                    std::cout << "Node " << nodeId << ": ";
                    for (const auto& block : blocks) {
                        std::cout << "[" << block[0] << "," << block[1] << "," << block[2] << "] ";
                    }
                    std::cout << "\n";
                }
                                
                std::cout << "Original: ...\n" 
                         << strOriginal << "\n";
                std::cout << "Fixed:    ...\n" 
                         << strFixed << "\n";
                
                // Print internal representations
                auto vecOrig = msaOriginal.getMSAVec();
                auto vecFixed = msaFixed.getMSAVec();
                
                std::cout << "\n--- Internal Representations ---\n";
                for (size_t i = 0; i < saveList.size(); i++) {
                    if (!saveList[i]) continue;
                    std::cout << "Node " << i << ":\n";
                    std::cout << "  Original: ";
                    for (auto v : vecOrig[i]) std::cout << v << " ";
                    std::cout << "\n  Fixed:    ";
                    for (auto v : vecFixed[i]) std::cout << v << " ";
                    std::cout << "\n";
                }
                
                return 1; // Stop at first mismatch for debugging
            } else {
                std::cout << "✓ Match (Orig len: " << msaOriginal.getMSAlength() 
                         << ", Fixed len: " << msaFixed.getMSAlength() << ")\n";
            }
            
        } catch (const std::exception& e) {
            std::cout << "✗ Exception: " << e.what() << "\n";
            return 1;
        }
    }
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "FINAL RESULTS: " << (numTrials - numMismatches) << "/" << numTrials 
              << " trials matched\n";
    
    if (numMismatches == 0) {
        std::cout << "✓ ALL TRIALS MATCHED!\n";
        return 0;
    } else {
        std::cout << "✗ " << numMismatches << " MISMATCHES FOUND\n";
        return 1;
    }
}