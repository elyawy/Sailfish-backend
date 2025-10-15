#include "../../../src/Simulator.h"
#include "../../../src/MSA.h"
#include "../../../src/MsaFixed.h"

int main() {
    // Simple tree with 3 leaves
    tree tree_("(A:0.1,B:0.2,C:0.3);", false);
    
    vector<DiscreteDistribution*> insertionDists(tree_.getNodesNum() - 1);
    vector<DiscreteDistribution*> deletionDists(tree_.getNodesNum() - 1);

    // Simple Zipf distribution
    DiscreteDistribution d1({0.5, 0.3, 0.2});
    
    fill(insertionDists.begin(), insertionDists.end(), &d1);
    fill(deletionDists.begin(), deletionDists.end(), &d1);

    vector<double> insertionRates(tree_.getNodesNum() - 1);
    vector<double> deletionRates(tree_.getNodesNum() - 1);

    fill(insertionRates.begin(), insertionRates.end(), 1.0);
    fill(deletionRates.begin(), deletionRates.end(), 1.0);

    SimulationProtocol protocol(&tree_);
    protocol.setInsertionLengthDistributions(insertionDists);
    protocol.setDeletionLengthDistributions(deletionDists);
    protocol.setInsertionRates(insertionRates);
    protocol.setDeletionRates(deletionRates);
    protocol.setSequenceSize(5);
    protocol.setSeed(42);

    Simulator sim(&protocol);
    auto saveList = sim.getNodesSaveList();

    // Generate one simulation
    auto blockmap = sim.generateSimulation();

    // Build with original implementation
    MSA msaOriginal(blockmap, tree_.getRoot(), saveList);
    std::string strOriginal = msaOriginal.generateMsaString();
    std::cout << "MSaOriginalLength=" << msaOriginal.getMSAlength() << "\n";
    // Build with FixedList implementation
    MsaFixed msaFixed(blockmap, tree_.getRoot(), saveList);

    std::string strFixed = msaFixed.generateMsaString();

    // Compare
    if (strOriginal == strFixed) {
        std::cout << "✓ MATCH! Both implementations produce identical MSA\n";
        std::cout << "MSA length: " << msaOriginal.getMSAlength() << "\n";
        return 0;
    } else {
        std::cout << "✗ MISMATCH! Implementations differ\n";
        std::cout << "Original MSA length: " << msaOriginal.getMSAlength() << "\n";
        std::cout << "Fixed MSA length: " << msaFixed.getMSAlength() << "\n";
        std::cout << "\n--- Original MSA ---\n" << strOriginal << "\n";
        std::cout << "\n--- Fixed MSA ---\n" << strFixed << "\n";
        return 1;
    }
}