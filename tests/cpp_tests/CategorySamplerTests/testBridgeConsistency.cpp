#include <iostream>
#include <vector>
#include <random>
#include "../../../src/CategorySampler.h"

int main() {
    const std::vector<std::vector<double>> matrix = {{0.01, 0.99},
                                                     {0.99, 0.01}};
    const std::vector<double> stationary = {0.5, 0.5};
    const size_t maxPathLength = 10;
    const int chainLength = 10000000;

    CategorySampler sampler(matrix, stationary, maxPathLength);
    std::mt19937_64 rng(42);

    // Build a long Markov chain
    std::vector<int> chain(chainLength);
    chain[0] = 0;
    for (int i = 1; i < chainLength; ++i) {
        chain[i] = sampler.drawSample(rng, chain[i-1]);
    }

    auto testConsistency = [&](int left, int right, size_t length) {
        // Count middle state distributions from the chain
        std::vector<int> chainCounts(2, 0);
        int chainTotal = 0;
        for (int i = 0; i < chainLength - (int)length - 1; ++i) {
            if (chain[i] == left && chain[i + length + 1] == right) {
                chainCounts[chain[i + 1]]++;
                chainTotal++;
            }
        }

        // Count middle state distributions from bridge sampling
        const int nBridgeTrials = 100000;
        std::vector<int> bridgeCounts(2, 0);
        for (int i = 0; i < nBridgeTrials; ++i) {
            auto path = sampler.sampleBridge(left, right, length, rng);
            bridgeCounts[path[0]]++;
        }

        std::cout << "Bridge " << left << "->?->" << right 
                  << " (length=" << length << "):\n";
        if (chainTotal == 0) {
            std::cout << "  Chain: no occurrences found\n";
        } else {
            std::cout << "  Chain:  P(0)=" << (double)chainCounts[0]/chainTotal 
                      << " P(1)=" << (double)chainCounts[1]/chainTotal 
                      << " (n=" << chainTotal << ")\n";
        }
        std::cout << "  Bridge: P(0)=" << (double)bridgeCounts[0]/nBridgeTrials 
                  << " P(1)=" << (double)bridgeCounts[1]/nBridgeTrials << "\n";
    };

    std::cout << "=== Length 1 bridges ===\n";
    testConsistency(0, 0, 1);
    testConsistency(1, 1, 1);
    testConsistency(0, 1, 1);
    testConsistency(1, 0, 1);

    std::cout << "\n=== Length 2 bridges ===\n";
    testConsistency(0, 0, 2);
    testConsistency(0, 1, 2);

    return 0;
}
