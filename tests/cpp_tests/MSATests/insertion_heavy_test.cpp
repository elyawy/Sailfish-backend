#include "../../../src/Simulator.h"
#include "../../../src/MSA.h"
#include "../../../src/MsaFixed.h"

int main() {
    // Tree with 2 children from root - allows overlapping insertions
    tree tree_("(A:0.1,(B:0.1,C:0.1):0.1);", false);
    
    vector<DiscreteDistribution*> insertionDists(tree_.getNodesNum() - 1);
    vector<DiscreteDistribution*> deletionDists(tree_.getNodesNum() - 1);

    // Distribution that always gives length 1
    DiscreteDistribution d1({1.0});
    
    fill(insertionDists.begin(), insertionDists.end(), &d1);
    fill(deletionDists.begin(), deletionDists.end(), &d1);

    vector<double> insertionRates(tree_.getNodesNum() - 1);
    vector<double> deletionRates(tree_.getNodesNum() - 1);

    // High insertion rate to guarantee insertions, no deletions
    fill(insertionRates.begin(), insertionRates.end(), 5.0);
    fill(deletionRates.begin(), deletionRates.end(), 1.0);

    SimulationProtocol protocol(&tree_);
    protocol.setInsertionLengthDistributions(insertionDists);
    protocol.setDeletionLengthDistributions(deletionDists);
    protocol.setInsertionRates(insertionRates);
    protocol.setDeletionRates(deletionRates);
    protocol.setSequenceSize(10);  // Small root sequence

    Simulator sim(&protocol);
    auto saveList = sim.getNodesSaveList();

    // Run multiple simulations to catch different insertion scenarios
    for (int trial = 0; trial < 1000; trial++) {
        protocol.setSeed(trial);

        auto blockmap = sim.generateSimulation();
        
        // Print the blockmaps to see where insertions occur
        std::cout << "\n=== Trial " << trial << " ===\n";
        for (const auto& [nodeId, blockTuple] : blockmap) {
            const BlockList& blocks = std::get<static_cast<int>(BLOCKLIST::BLOCKS)>(blockTuple);
            std::cout << "Node " << nodeId << " blocks:\n";
            for (const auto& block : blocks) {
                std::cout << "  [pos=" << block[0] 
                         << ", len=" << block[1] 
                         << ", ins=" << block[2] << "]\n";
            }
        }

        // Build with both implementations
        MSA msaOriginal(blockmap, tree_.getRoot(), saveList);
        MsaFixed msaFixed(blockmap, tree_.getRoot(), saveList);

        std::string strOriginal = msaOriginal.generateMsaString();
        std::string strFixed = msaFixed.generateMsaString();

        // Compare
        if (strOriginal != strFixed) {
            std::cout << "\n✗ MISMATCH found!\n";
            std::cout << "Original MSA length: " << msaOriginal.getMSAlength() << "\n";
            std::cout << "Fixed MSA length: " << msaFixed.getMSAlength() << "\n";
            
            std::cout << "\n--- Original MSA ---\n" << strOriginal << "\n";
            std::cout << "\n--- Fixed MSA ---\n" << strFixed << "\n";
            
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
            
            return 1;  // Stop at first mismatch
        } else {
            std::cout << "✓ Match for trial " << trial << "\n";
        }
    }

    std::cout << "\n✓ ALL TRIALS MATCHED!\n";
    return 0;
}