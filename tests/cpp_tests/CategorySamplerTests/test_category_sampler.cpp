#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include "../../../src/CategorySampler.h"

void testIndependentSampling() {
    std::cout << "=== Test 1: Independent Sampling (correlation = 0) ===" << std::endl;
    
    // 4 categories with equal probabilities
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    CategorySampler sampler(probs, 0.0);
    
    // Sample 10000 sites
    std::map<int, int> counts;
    const int numSamples = 10000;
    
    for (int i = 0; i < numSamples; ++i) {
        int category = sampler.drawSample();
        counts[category]++;
    }
    
    std::cout << "Expected frequency: 0.25 for each category" << std::endl;
    std::cout << "Observed frequencies:" << std::endl;
    for (const auto& pair : counts) {
        double freq = static_cast<double>(pair.second) / numSamples;
        std::cout << "  Category " << pair.first << ": " << freq 
                  << " (count: " << pair.second << ")" << std::endl;
    }
    std::cout << std::endl;
}

void testPerfectCorrelation() {
    std::cout << "=== Test 2: Perfect Correlation (correlation = 1.0) ===" << std::endl;
    
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    CategorySampler sampler(probs, 1.0);
    
    // Sample 100 sites - should all be the same
    int firstCategory = sampler.drawSample();
    bool allSame = true;
    
    std::cout << "First category: " << firstCategory << std::endl;
    std::cout << "Next 99 categories: ";
    
    for (int i = 0; i < 99; ++i) {
        int category = sampler.drawSample();
        if (i < 10) std::cout << category << " ";
        if (category != firstCategory) {
            allSame = false;
        }
    }
    
    std::cout << std::endl;
    std::cout << "All categories same as first? " << (allSame ? "YES ✓" : "NO ✗") << std::endl;
    std::cout << std::endl;
}

void testModerateCorrelation() {
    std::cout << "=== Test 3: Moderate Correlation (correlation = 0.7) ===" << std::endl;
    
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    CategorySampler sampler(probs, 0.7);
    
    // Count transitions
    int numTransitions = 0;
    int previousCategory = sampler.drawSample();
    const int numSamples = 1000;
    
    for (int i = 0; i < numSamples - 1; ++i) {
        int currentCategory = sampler.drawSample();
        if (currentCategory != previousCategory) {
            numTransitions++;
        }
        previousCategory = currentCategory;
    }
    
    double transitionRate = static_cast<double>(numTransitions) / (numSamples - 1);
    double expectedTransitionRate = (1.0 - 0.7) * (1.0 - 0.25);  // 0.3 * 0.75 = 0.225
    
    
    std::cout << "Number of transitions: " << numTransitions << " out of " << (numSamples - 1) << std::endl;
    std::cout << "Observed transition rate: " << transitionRate << std::endl;
    std::cout << "Expected transition rate: ~" << expectedTransitionRate << std::endl;

    bool withinRange = std::abs(transitionRate - expectedTransitionRate) < 0.02;
    std::cout << "Within expected range? " << (withinRange ? "YES ✓" : "NO ✗") << std::endl;
    std::cout << std::endl;
}

void testNonUniformDistribution() {
    std::cout << "=== Test 4: Non-uniform Distribution ===" << std::endl;
    
    // Heavily skewed distribution
    std::vector<MDOUBLE> probs = {0.1, 0.2, 0.3, 0.4};
    CategorySampler sampler(probs, 0.0);  // Independent sampling
    
    std::map<int, int> counts;
    const int numSamples = 10000;
    
    for (int i = 0; i < numSamples; ++i) {
        int category = sampler.drawSample();
        counts[category]++;
    }
    
    std::cout << "Expected vs Observed frequencies:" << std::endl;
    for (size_t i = 0; i < probs.size(); ++i) {
        double observed = static_cast<double>(counts[i]) / numSamples;
        std::cout << "  Category " << i << ": expected=" << probs[i] 
                  << ", observed=" << observed << std::endl;
    }
    std::cout << std::endl;
}

void testReset() {
    std::cout << "=== Test 5: Reset Functionality ===" << std::endl;
    
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    CategorySampler sampler(probs, 1.0);  // Perfect correlation
    
    int firstSeq = sampler.drawSample();
    std::cout << "First sequence starts with category: " << firstSeq << std::endl;
    
    for (int i = 0; i < 10; ++i) {
        sampler.drawSample();  // Should all be same
    }
    
    sampler.reset();
    int secondSeq = sampler.drawSample();
    std::cout << "After reset, new sequence starts with category: " << secondSeq << std::endl;
    std::cout << "Categories differ after reset? " << (firstSeq != secondSeq ? "Possibly (good)" : "No (might be same by chance)") << std::endl;
    std::cout << std::endl;
}

void testWithInvariantSites() {
    std::cout << "=== Test 6: With Invariant Sites ===" << std::endl;
    
    // 4 gamma categories + 1 invariant category
    double invariantProp = 0.2;
    std::vector<MDOUBLE> probs = {
        0.25 * (1 - invariantProp),  // 0.2
        0.25 * (1 - invariantProp),  // 0.2
        0.25 * (1 - invariantProp),  // 0.2
        0.25 * (1 - invariantProp),  // 0.2
        invariantProp               // 0.20 (invariant)
    };
    
    CategorySampler sampler(probs, 0.0);
    
    std::map<int, int> counts;
    const int numSamples = 10000;
    
    for (int i = 0; i < numSamples; ++i) {
        int category = sampler.drawSample();
        counts[category]++;
    }
    
    std::cout << "Expected vs Observed frequencies (with invariant):" << std::endl;
    for (size_t i = 0; i < probs.size(); ++i) {
        double observed = static_cast<double>(counts[i]) / numSamples;
        std::string label = (i == 4) ? " (invariant)" : "";
        std::cout << "  Category " << i << label << ": expected=" << probs[i] 
                  << ", observed=" << observed << std::endl;
    }
    std::cout << std::endl;
}

void testNonUniformModerateCorrelation() {
    std::cout << "=== Test 7: Non-uniform Distribution with Moderate Correlation ===" << std::endl;
    
    // Heavily skewed distribution
    std::vector<MDOUBLE> probs = {0.1, 0.2, 0.3, 0.4};
    double correlation = 0.6;
    CategorySampler sampler(probs, correlation);
    
    const int numSamples = 10000000;  // Larger sample for better statistics
    
    // Count category frequencies
    std::map<int, int> counts;
    // Count transitions from each category
    std::map<int, int> fromCounts;
    std::map<int, std::map<int, int>> transitionCounts;
    
    int previousCategory = sampler.drawSample();
    counts[previousCategory]++;
    // fromCounts[previousCategory]++;
    
    for (int i = 0; i < numSamples - 1; ++i) {
        int currentCategory = sampler.drawSample();
        counts[currentCategory]++;
        transitionCounts[previousCategory][currentCategory]++;
        fromCounts[previousCategory]++;
        previousCategory = currentCategory;
    }
    
    std::cout << "\n1. Stationary Distribution Check:" << std::endl;
    std::cout << "   (Long-run frequencies should match π)" << std::endl;
    for (size_t i = 0; i < probs.size(); ++i) {
        double observed = static_cast<double>(counts[i]) / numSamples;
        std::cout << "   Category " << i << ": expected=" << probs[i] 
                  << ", observed=" << observed;
        bool match = std::abs(observed - probs[i]) < 0.01;
        std::cout << (match ? " ✓" : " ✗") << std::endl;
    }
    
    std::cout << "\n2. Transition Probabilities:" << std::endl;
    std::cout << "   P[i][j] = ρ * δ(i,j) + (1-ρ) * π[j]" << std::endl;
    std::cout << "   where ρ = " << correlation << std::endl;
    
    for (size_t i = 0; i < probs.size(); ++i) {
        if (fromCounts[i] == 0) continue;
        
        std::cout << "\n   From Category " << i << ":" << std::endl;
        for (size_t j = 0; j < probs.size(); ++j) {
            double observedProb = static_cast<double>(transitionCounts[i][j]) / fromCounts[i];
            double expectedProb;
            if (i == j) {
                expectedProb = correlation + (1.0 - correlation) * probs[j];
            } else {
                expectedProb = (1.0 - correlation) * probs[j];
            }
            
            std::cout << "     → Category " << j << ": ";
            std::cout << "expected=" << expectedProb << ", observed=" << observedProb;
            
            bool match = std::abs(observedProb - expectedProb) < 0.01;
            std::cout << (match ? " ✓" : " ✗") << std::endl;
        }
    }
    
    std::cout << "\n3. Overall Transition Rate:" << std::endl;
    int totalTransitions = 0;
    for (size_t i = 0; i < probs.size(); ++i) {
        for (size_t j = 0; j < probs.size(); ++j) {
            if (i != j) {
                totalTransitions += transitionCounts[i][j];
            }
        }
    }
    
    double observedTransitionRate = static_cast<double>(totalTransitions) / (numSamples - 1);
    
    // Expected transition rate for non-uniform distribution
    // E[P(switch)] = Σ_i π[i] * (1 - P[i][i])
    //              = Σ_i π[i] * (1 - ρ - (1-ρ)*π[i])
    //              = Σ_i π[i] * (1-ρ) * (1-π[i])
    //              = (1-ρ) * Σ_i π[i] * (1-π[i])
    //              = (1-ρ) * (1 - Σ_i π[i]²)
    double sumPiSquared = 0.0;
    for (size_t i = 0; i < probs.size(); ++i) {
        sumPiSquared += probs[i] * probs[i];
    }
    double expectedTransitionRate = (1.0 - correlation) * (1.0 - sumPiSquared);
    
    std::cout << "   Observed: " << observedTransitionRate << std::endl;
    std::cout << "   Expected: " << expectedTransitionRate << std::endl;
    bool match = std::abs(observedTransitionRate - expectedTransitionRate) < 0.01;
    std::cout << "   Match? " << (match ? "YES ✓" : "NO ✗") << std::endl;
    
    std::cout << "\n4. Auto-correlation Check:" << std::endl;
    std::cout << "   More probable categories should have longer runs" << std::endl;
    
    for (size_t i = 0; i < probs.size(); ++i) {
        if (fromCounts[i] == 0) continue;
        
        // P(stay in i | currently in i)
        double stayProb = static_cast<double>(transitionCounts[i][i]) / fromCounts[i];
        double expectedStayProb = correlation + (1.0 - correlation) * probs[i];
        
        // Expected run length = 1 / P(leave)
        double expectedRunLength = 1.0 / (1.0 - expectedStayProb);
        double observedRunLength = 1.0 / (1.0 - stayProb);
        
        std::cout << "   Category " << i << " (π=" << probs[i] << "): ";
        std::cout << "expected run length=" << expectedRunLength 
                  << ", observed=" << observedRunLength << std::endl;
    }
    
    std::cout << std::endl;
}


void testAliasMethodAccuracy() {
    std::cout << "=== Alias Method Accuracy Test ===" << std::endl;
    
    std::vector<double> probs = {0.25, 0.25, 0.25, 0.25};
    DiscreteDistribution sampler(probs);
    
    std::map<int, int> counts;
    const int numSamples = 1000000;  // 1 million samples
    
    for (int i = 0; i < numSamples; ++i) {
        int sample = sampler.drawSample();
        counts[sample]++;
    }
    
    std::cout << "After " << numSamples << " samples:" << std::endl;
    for (const auto& pair : counts) {
        double freq = static_cast<double>(pair.second) / numSamples;
        double expected = 0.25;
        double error = std::abs(freq - expected);
        std::cout << "  Value " << pair.first << ": " << freq 
                  << " (error: " << error << ")" << std::endl;
    }
}

int main() {
    std::cout << "CategorySampler Test Suite" << std::endl;
    std::cout << "===========================" << std::endl << std::endl;
    
    try {
        testIndependentSampling();
        testPerfectCorrelation();
        testModerateCorrelation();
        testNonUniformDistribution();
        testReset();
        testWithInvariantSites();
        testNonUniformModerateCorrelation();
        testAliasMethodAccuracy();


        std::cout << "All tests completed successfully! ✓" << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}