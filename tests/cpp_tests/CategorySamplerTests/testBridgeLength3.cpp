#include <iostream>
#include <vector>
#include <random>
#include "../../../src/CategorySampler.h"

int main() {
    const std::vector<std::vector<double>> matrix = {{0.01, 0.99},
                                                     {0.99, 0.01}};
    const std::vector<double> stationary = {0.5, 0.5};
    const size_t maxPathLength = 10;
    const int chainLength = 100000000; // 100M for rare combinations

    CategorySampler sampler(matrix, stationary, maxPathLength);
    std::mt19937_64 rng(42);

    // Build a long Markov chain
    std::cout << "Building chain of length " << chainLength << "...\n";
    std::vector<int> chain(chainLength);
    chain[0] = 0;
    for (int i = 1; i < chainLength; ++i) {
        chain[i] = sampler.drawSample(rng, chain[i-1]);
    }
    std::cout << "Done.\n\n";

    auto testConsistency = [&](int left, int right, size_t length) {
        int nStates = 2;
        // For each intermediate position, count state frequencies
        std::vector<std::vector<int>> chainCounts(length, std::vector<int>(nStates, 0));
        int chainTotal = 0;
        for (int i = 0; i < chainLength - (int)length - 1; ++i) {
            if (chain[i] == left && chain[i + length + 1] == right) {
                for (size_t pos = 0; pos < length; ++pos) {
                    chainCounts[pos][chain[i + 1 + pos]]++;
                }
                chainTotal++;
            }
        }

        // Bridge sampling
        const int nBridgeTrials = 100000;
        std::vector<std::vector<int>> bridgeCounts(length, std::vector<int>(nStates, 0));
        for (int i = 0; i < nBridgeTrials; ++i) {
            auto path = sampler.sampleBridge(left, right, length, rng);
            for (size_t pos = 0; pos < length; ++pos) {
                bridgeCounts[pos][path[pos]]++;
            }
        }

        std::cout << "Bridge " << left << "->???->" << right 
                  << " (length=" << length << "):\n";
        if (chainTotal == 0) {
            std::cout << "  Chain: no occurrences found\n";
        } else {
            std::cout << "  Chain (n=" << chainTotal << "):\n";
            for (size_t pos = 0; pos < length; ++pos) {
                std::cout << "    pos " << pos+1 << ": P(0)=" << (double)chainCounts[pos][0]/chainTotal
                          << " P(1)=" << (double)chainCounts[pos][1]/chainTotal << "\n";
            }
        }
        std::cout << "  Bridge:\n";
        for (size_t pos = 0; pos < length; ++pos) {
            std::cout << "    pos " << pos+1 << ": P(0)=" << (double)bridgeCounts[pos][0]/nBridgeTrials
                      << " P(1)=" << (double)bridgeCounts[pos][1]/nBridgeTrials << "\n";
        }
    };

    std::cout << "=== Length 3 bridges ===\n";
    testConsistency(0, 0, 3);
    testConsistency(1, 1, 3);
    testConsistency(0, 1, 3);
    testConsistency(1, 0, 3);

    return 0;
}
