#include <iostream>
#include <random>
#include <cassert>

#include "../../../src/CategorySampler.h"

int main() {
    // Transition matrix
    std::vector<std::vector<double>> transitionMatrix = {
        {5.87080948e-01, 2.38313141e-01, 1.05565146e-01, 4.51370352e-02, 1.73360988e-02, 5.40148139e-03, 1.09522164e-03, 7.09272985e-05},
        {2.38313141e-01, 2.88960901e-01, 2.15507513e-01, 1.37139486e-01, 7.52696478e-02, 3.35508994e-02, 1.01631898e-02, 1.09522164e-03},
        {1.05565146e-01, 2.15507513e-01, 2.27831989e-01, 1.92946850e-01, 1.38204384e-01, 8.09917374e-02, 3.35508994e-02, 5.40148139e-03},
        {4.51370352e-02, 1.37139486e-01, 1.92946850e-01, 2.07576760e-01, 1.86389738e-01, 1.38204384e-01, 7.52696478e-02, 1.73360988e-02},
        {1.73360988e-02, 7.52696478e-02, 1.38204384e-01, 1.86389738e-01, 2.07576760e-01, 1.92946850e-01, 1.37139486e-01, 4.51370352e-02},
        {5.40148139e-03, 3.35508994e-02, 8.09917374e-02, 1.38204384e-01, 1.92946850e-01, 2.27831989e-01, 2.15507513e-01, 1.05565146e-01},
        {1.09522164e-03, 1.01631898e-02, 3.35508994e-02, 7.52696478e-02, 1.37139486e-01, 2.15507513e-01, 2.88960901e-01, 2.38313141e-01},
        {7.09272985e-05, 1.09522164e-03, 5.40148139e-03, 1.73360988e-02, 4.51370352e-02, 1.05565146e-01, 2.38313141e-01, 5.87080948e-01}
    };
    
    // Stationary probabilities (uniform for simplicity)
    std::vector<double> stationaryProbs(8, 1.0/8.0);
    
    // Create sampler with max path length
    size_t maxPathLength = 100;
    CategorySampler sampler(transitionMatrix, stationaryProbs, maxPathLength);
    
    std::mt19937_64 rng(42);
    
    std::cout << "=== Test 1: Bridge of length 50 ===" << std::endl;
    auto bridge1 = sampler.sampleBridge(0, 7, 50, rng);
    assert(bridge1.size() == 50);
    std::cout << "Bridge from state 0 to state 7, length 50:" << std::endl;
    std::cout << "First 10: ";
    for (size_t i = 0; i < 10; ++i) std::cout << bridge1[i] << " ";
    std::cout << std::endl;
    std::cout << "Last 10: ";
    for (size_t i = 40; i < 50; ++i) std::cout << bridge1[i] << " ";
    std::cout << std::endl << std::endl;
    
    std::cout << "=== Test 2: Nested bridge of length 10 ===" << std::endl;
    size_t leftIdx = 20;
    size_t rightIdx = 30;
    size_t leftState = bridge1[leftIdx];
    size_t rightState = bridge1[rightIdx];
    auto bridge2 = sampler.sampleBridge(leftState, rightState, 10, rng);
    assert(bridge2.size() == 10);
    std::cout << "Nested bridge from state " << leftState << " to state " << rightState << ", length 10:" << std::endl;
    for (size_t i = 0; i < 10; ++i) std::cout << bridge2[i] << " ";
    std::cout << std::endl << std::endl;
    
    std::cout << "=== Test 3: Edge case - length 1 bridge ===" << std::endl;
    auto bridge3 = sampler.sampleBridge(2, 5, 1, rng);
    assert(bridge3.size() == 1);
    std::cout << "Bridge from state 2 to state 5, length 1: " << bridge3[0] << std::endl << std::endl;
    
    std::cout << "=== Test 4: Same start and end - length 1 ===" << std::endl;
    auto bridge4 = sampler.sampleBridge(3, 3, 1, rng);
    assert(bridge4.size() == 1);
    std::cout << "Bridge from state 3 to state 3, length 1: " << bridge4[0] << std::endl << std::endl;
    
    std::cout << "=== Test 5: Maximum length bridge ===" << std::endl;
    auto bridge5 = sampler.sampleBridge(1, 6, maxPathLength, rng);
    assert(bridge5.size() == maxPathLength);
    std::cout << "Bridge from state 1 to state 6, length " << maxPathLength << ":" << std::endl;
    std::cout << "First 10: ";
    for (size_t i = 0; i < 10; ++i) std::cout << bridge5[i] << " ";
    std::cout << std::endl;
    std::cout << "Last 10: ";
    for (size_t i = maxPathLength - 10; i < maxPathLength; ++i) std::cout << bridge5[i] << " ";
    std::cout << std::endl << std::endl;
    
    std::cout << "All tests passed!" << std::endl;
    
    return 0;
}