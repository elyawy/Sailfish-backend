#include <iostream>
#include <vector>
#include <random>
#include "../../../src/CategorySampler.h"

int main() {
    const std::vector<std::vector<double>> matrix = {{0.01, 0.99},
                                                     {0.99, 0.01}};
    const std::vector<double> stationary = {0.5, 0.5};
    const size_t maxPathLength = 100;
    const int nTrials = 100000;

    CategorySampler sampler(matrix, stationary, maxPathLength);
    std::mt19937_64 rng(42);

    auto runTest = [&](size_t left, size_t right, size_t length) {
        std::vector<int> counts(2, 0);
        for (int i = 0; i < nTrials; ++i) {
            auto path = sampler.sampleBridge(left, right, length, rng);
            counts[path[0]]++;
        }
        std::cout << "Bridge " << left << "->?->" << right 
                  << " (length=" << length << "): "
                  << "P(0)=" << (double)counts[0]/nTrials 
                  << " P(1)=" << (double)counts[1]/nTrials << "\n";
    };

    std::cout << "=== Length 1 bridges ===\n";
    runTest(0, 0, 1);   // must be 1
    runTest(1, 1, 1);   // must be 0
    runTest(0, 1, 1);   // ~50/50
    runTest(1, 0, 1);   // ~50/50

    std::cout << "\n=== Length 2 bridges ===\n";
    runTest(0, 0, 2);   // first state: ~50/50 conditioned on returning to 0
    runTest(0, 1, 2);   // first state: strongly biased toward 1

    return 0;
}
