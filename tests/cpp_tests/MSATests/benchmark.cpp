#include <chrono>
#include <iomanip>
#include "../../../src/Simulator.h"
#include "../../../src/MSA.h"
#include "../../../src/MsaFixed.h"

int main() {
    // Fixed parameters
    std::string newick = "../../trees/normalbranches_nLeaves100.treefile";
    int rootLength = 100;
    double insertionRate = 0.05;
    double deletionRate = 0.05;
    int numTrials = 1000;
    
    std::cout << "=== MSA Implementation Benchmark ===\n";
    std::cout << "Tree: " << newick << "\n";
    std::cout << "Root length: " << rootLength << "\n";
    std::cout << "Insertion rate: " << insertionRate << "\n";
    std::cout << "Deletion rate: " << deletionRate << "\n";
    std::cout << "Number of trials: " << numTrials << "\n\n";
    
    tree tree_(newick, true);
    
    vector<DiscreteDistribution*> insertionDists(tree_.getNodesNum() - 1);
    vector<DiscreteDistribution*> deletionDists(tree_.getNodesNum() - 1);
    
    // Zipf distribution with 'a' parameter of 1.7
    DiscreteDistribution indelDist({0.500597013622622,0.15407680416466002,0.07733607404815325,0.04742269916833997,0.03245183409648796,0.023802968878603163,0.018315618393686197,0.014596047786709592,0.011947471092366838,0.009988223562671806,0.008494164935846923,0.007326223038978798,0.00639417555007358,0.005637292815586038,0.005013408742738382,0.004492460672380736,0.004052519265065782,0.0036772655722418293,0.0033543403993036617,0.0030742364096058624,0.00282953749580224,0.002614385926699555,0.0024241022418638207,0.0022549096413393907,0.0021037311601300934,0.001968038376604936,0.0018457371577210268,0.001735080389074546,0.001634600615111441,0.00154305754135923,0.001459396750517745,0.0013827169647433223,0.0013122438819642994,0.0012473091132895959,0.0011873331108893934,0.0011318112414127814,0.0010803023567287186,0.0010324193607649647,0.0009878213819912308,0.000946207245246948,0.0009073100010352715,0.000870892320076723,0.0008367425994772478,0.0008046716569959321,0.0007745099135856162,0.0007461049831118826,0.0007193196030519133,0.000694029851883145,0.0006701236084401137,0.000647499216243446,0.0006260643220715613,0.0006057348631513914,0.0005864341815220788,0.0005680922475578417,0.0005506449774670669,0.0005340336319283628,0.0005182042849716617,0.0005031073538361537,0.0004886971818952107,0.00047493166787853947,0.0004617719355816668,0.00044918203906341153,0.00043712869901847034,0.00042558106659524523,0.0004145105114255328,0.0004038904310565982,0.0003936960793390335,0.00038390441163517177,0.00037449394498066375,0.0003654446315627292,0.0003567377440781168,0.00034835577170658,0.00034028232558560196,0.0003325020528024772,0.00032500055803344507,0.000317764332058729,0.0003107806864690648,0.0003040376939552896,0.00029752413363926596,0.00029122944096306576,0.0002851436617049981,0.00027925740973663626,0.00027356182817527154,0.0002680485536218626,0.00026270968320613465,0.0002575377441885291,0.00025252566589363025,0.00024766675377188466,0.00024295466540621684,0.00023838338829779629,0.00023394721928100186,0.00022964074543174937,0.0002254588263460113,0.00022139657767671146,0.00021744935582738057,0.00021361274371013342,0.00020988253748378816,0.00020625473419539554,0.00020272552025516334,0.0001992912606808326});
    
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
    
    Simulator sim(&protocol);
    auto saveList = sim.getNodesSaveList();
    
    // Pre-generate all blockmaps to avoid including simulation time
    std::cout << "Pre-generating " << numTrials << " blockmaps...\n";
    std::vector<BlockMap> blockmaps;
    blockmaps.reserve(numTrials);
    
    for (int trial = 0; trial < numTrials; trial++) {
        protocol.setSeed(trial);
        blockmaps.push_back(sim.generateSimulation());
    }
    std::cout << "Done!\n\n";
    
    // Benchmark Original MSA (SuperSequence) - construction only
    std::cout << "Benchmarking Original MSA (SuperSequence) construction...\n";
    auto startOriginal = std::chrono::high_resolution_clock::now();
    
    for (int trial = 0; trial < numTrials; trial++) {
        MSA msaOriginal(blockmaps[trial], tree_.getRoot(), saveList);
    }
    
    auto endOriginal = std::chrono::high_resolution_clock::now();
    auto durationOriginal = std::chrono::duration_cast<std::chrono::milliseconds>(endOriginal - startOriginal);
    
    std::cout << "Original MSA total time: " << durationOriginal.count() << " ms\n";
    std::cout << "Original MSA avg per trial: " 
              << std::fixed << std::setprecision(3) 
              << (durationOriginal.count() / (double)numTrials) << " ms\n\n";
    
    // Benchmark Fixed MSA (FixedList) - construction only
    std::cout << "Benchmarking Fixed MSA (FixedList) construction...\n";
    auto startFixed = std::chrono::high_resolution_clock::now();
    
    for (int trial = 0; trial < numTrials; trial++) {
        MsaFixed msaFixed(blockmaps[trial], tree_.getRoot(), saveList);
    }
    
    auto endFixed = std::chrono::high_resolution_clock::now();
    auto durationFixed = std::chrono::duration_cast<std::chrono::milliseconds>(endFixed - startFixed);
    
    std::cout << "Fixed MSA total time: " << durationFixed.count() << " ms\n";
    std::cout << "Fixed MSA avg per trial: " 
              << std::fixed << std::setprecision(3) 
              << (durationFixed.count() / (double)numTrials) << " ms\n\n";
    
    // Calculate speedup
    double speedup = (double)durationOriginal.count() / (double)durationFixed.count();
    
    std::cout << std::string(60, '=') << "\n";
    std::cout << "RESULTS:\n";
    std::cout << "Original MSA: " << durationOriginal.count() << " ms\n";
    std::cout << "Fixed MSA:    " << durationFixed.count() << " ms\n";
    std::cout << "Speedup:      " << std::fixed << std::setprecision(2) << speedup << "x ";
    
    if (speedup > 1.0) {
        std::cout << "(Fixed is faster)\n";
    } else if (speedup < 1.0) {
        std::cout << "(Original is faster)\n";
    } else {
        std::cout << "(Same speed)\n";
    }
    std::cout << std::string(60, '=') << "\n";
    
    return 0;
}