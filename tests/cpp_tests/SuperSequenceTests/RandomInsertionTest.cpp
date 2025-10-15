#include <iostream>
#include <vector>
#include <cassert>
#include <sstream>
#include <random>
#include <chrono>
#include "../../src/FixedList.h"
#include "../../src/SuperSequence.h"

// Helper function to advance FixedList iterator manually
FixedList::iterator advanceFixed(FixedList::iterator it, int n) {
    for (int i = 0; i < n; ++i) {
        ++it;
    }
    return it;
}

// Helper function to advance SuperSequence iterator manually
SuperSequence::SequenceType::iterator advanceSuper(SuperSequence::SequenceType::iterator it, int n) {
    for (int i = 0; i < n; ++i) {
        ++it;
    }
    return it;
}

// Helper to get traversal as string for FixedList (skipping first element at index 0)
std::string getTraversalFixed(FixedList& fl) {
    std::stringstream ss;
    bool first = true;
    auto it = fl.begin();
    ++it;  // Skip element at index 0 to match SuperSequence
    
    for (; it != fl.end(); ++it) {
        if (!first) ss << "-";
        ss << *it;
        first = false;
    }
    return ss.str();
}

// Helper to get traversal as string for SuperSequence (position values)
std::string getTraversalSuper(SuperSequence& ss) {
    std::stringstream result;
    bool first = true;
    for (auto it = ss.begin(); it != ss.end(); ++it) {
        if (!first) result << "-";
        result << it->position;
        first = false;
    }
    return result.str();
}

void assertTraversalsMatch(FixedList& fl, SuperSequence& ss, const std::string& testName) {
    std::string fixedTraversal = getTraversalFixed(fl);
    std::string superTraversal = getTraversalSuper(ss);
    
    if (fixedTraversal != superTraversal) {
        std::cout << "FAILED: " << testName << "\n";
        std::cout << "  FixedList:     " << fixedTraversal << "\n";
        std::cout << "  SuperSequence: " << superTraversal << "\n";
        assert(false);
    }
}

void testRandomInsertions(int numInsertions, int initialSize, unsigned int seed) {
    std::cout << "\n=== Test: " << numInsertions << " random insertions (seed=" << seed << ") ===\n";
    
    // Calculate required capacity
    size_t requiredCapacity = 1 + initialSize + numInsertions;  // phantom + initial + insertions
    
    FixedList fl(requiredCapacity);
    fl.initialize(initialSize);
    
    SuperSequence ss(initialSize, 1);
    
    std::mt19937 rng(seed);
    
    std::cout << "Initial size: " << initialSize << "\n";
    std::cout << "Performing " << numInsertions << " random insertions...\n";
    
    assertTraversalsMatch(fl, ss, "Initial state");
    
    for (int i = 0; i < numInsertions; ++i) {
        // Current real sequence size (excluding phantom in FixedList)
        size_t currentSize = ss.size();
        
        // Generate random position to insert after (0 to currentSize-1)
        std::uniform_int_distribution<size_t> dist(0, currentSize - 1);
        size_t insertAfterPos = dist(rng);
        
        // Randomly decide if this should be a reference or not
        std::uniform_int_distribution<int> refDist(0, 9);
        bool isReference = (refDist(rng) == 0);  // 10% chance of reference
        
        // Insert in FixedList
        auto flIt = fl.begin();
        ++flIt;  // Skip phantom
        flIt = advanceFixed(flIt, insertAfterPos);
        fl.insertAfter(flIt, isReference);
        
        // Insert in SuperSequence
        auto ssIt = ss.begin();
        ssIt = advanceSuper(ssIt, insertAfterPos);
        ++ssIt;  // Move past the position to insert before
        size_t randomPos = (i == 0) ? ss.getRandomSequencePosition() : ss.incrementRandomSequencePosition();
        ss.insertItemAtPosition(ssIt, randomPos, isReference);
        
        // Verify every 1000 insertions
        if ((i + 1) % 1000 == 0) {
            std::stringstream name;
            name << "After " << (i + 1) << " insertions";
            assertTraversalsMatch(fl, ss, name.str());
            std::cout << "  ✓ Verified at " << (i + 1) << " insertions (size: " << ss.size() << ")\n";
        }
    }
    
    // Final verification
    assertTraversalsMatch(fl, ss, "Final state after all insertions");
    
    std::cout << "  Final FixedList size: " << fl.size() << " (including phantom)\n";
    std::cout << "  Final SuperSequence size: " << ss.size() << "\n";
    std::cout << "  ✓ ALL " << numInsertions << " INSERTIONS VERIFIED!\n";
}

void testRandomInsertionsWithReferences(int numInsertions, int initialSize, unsigned int seed) {
    std::cout << "\n=== Test: " << numInsertions << " random insertions with references (seed=" << seed << ") ===\n";
    
    size_t requiredCapacity = 1 + initialSize + numInsertions;
    
    FixedList fl(requiredCapacity);
    fl.initialize(initialSize);
    
    SuperSequence ss(initialSize, 1);
    
    std::mt19937 rng(seed);
    
    std::cout << "Initial size: " << initialSize << "\n";
    std::cout << "Performing " << numInsertions << " random insertions with references...\n";
    
    // Track which positions we've referenced
    std::vector<size_t> referencedPositions;
    
    for (int i = 0; i < numInsertions; ++i) {
        size_t currentSize = ss.size();
        
        std::uniform_int_distribution<size_t> dist(0, currentSize - 1);
        size_t insertAfterPos = dist(rng);
        
        // 30% chance of reference
        std::uniform_int_distribution<int> refDist(0, 9);
        bool isReference = (refDist(rng) < 3);
        
        // Insert in FixedList
        auto flIt = fl.begin();
        ++flIt;  // Skip phantom
        flIt = advanceFixed(flIt, insertAfterPos);
        size_t newIndexFL = *fl.insertAfter(flIt, isReference);
        
        // Insert in SuperSequence
        auto ssIt = ss.begin();
        ssIt = advanceSuper(ssIt, insertAfterPos);
        ++ssIt;
        size_t randomPos = (i == 0) ? ss.getRandomSequencePosition() : ss.incrementRandomSequencePosition();
        auto newSSIt = ss.insertItemAtPosition(ssIt, randomPos, isReference);
        
        if (isReference) {
            referencedPositions.push_back(i);
        }
        
        // Also randomly reference some existing positions
        if (currentSize > 5 && refDist(rng) < 2) {  // 20% chance
            std::uniform_int_distribution<size_t> existingDist(0, currentSize - 1);
            size_t existingPos = existingDist(rng);
            
            flIt = fl.begin();
            ++flIt;  // Skip phantom
            flIt = advanceFixed(flIt, existingPos);
            fl.referencePosition(flIt);
            
            ssIt = ss.begin();
            ssIt = advanceSuper(ssIt, existingPos);
            ss.referencePosition(ssIt);
        }
        
        if ((i + 1) % 1000 == 0) {
            std::stringstream name;
            name << "After " << (i + 1) << " insertions";
            assertTraversalsMatch(fl, ss, name.str());
            std::cout << "  ✓ Verified at " << (i + 1) << " insertions\n";
        }
    }
    
    // Final verification
    assertTraversalsMatch(fl, ss, "Final state after all insertions");
    
    // Set absolute positions and verify MSA lengths
    fl.setAbsolutePositions();
    ss.setAbsolutePositions();
    
    assert(fl.getMsaSequenceLength() == ss.getMsaSequenceLength());
    
    std::cout << "  Referenced positions: " << referencedPositions.size() << "\n";
    std::cout << "  Final MSA length: " << fl.getMsaSequenceLength() << "\n";
    std::cout << "  ✓ ALL " << numInsertions << " INSERTIONS WITH REFERENCES VERIFIED!\n";
}

void testStressTest(int numInsertions, unsigned int seed) {
    std::cout << "\n=== Stress Test: " << numInsertions << " insertions with small initial size (seed=" << seed << ") ===\n";
    
    int initialSize = 3;  // Start very small
    size_t requiredCapacity = 1 + initialSize + numInsertions;
    
    FixedList fl(requiredCapacity);
    fl.initialize(initialSize);
    
    SuperSequence ss(initialSize, 1);
    
    std::mt19937 rng(seed);
    
    std::cout << "Initial size: " << initialSize << "\n";
    std::cout << "Target insertions: " << numInsertions << "\n";
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < numInsertions; ++i) {
        size_t currentSize = ss.size();
        
        std::uniform_int_distribution<size_t> dist(0, currentSize - 1);
        size_t insertAfterPos = dist(rng);
        
        // Insert in FixedList
        auto flIt = fl.begin();
        ++flIt;
        flIt = advanceFixed(flIt, insertAfterPos);
        fl.insertAfter(flIt, false);
        
        // Insert in SuperSequence
        auto ssIt = ss.begin();
        ssIt = advanceSuper(ssIt, insertAfterPos);
        ++ssIt;
        size_t randomPos = (i == 0) ? ss.getRandomSequencePosition() : ss.incrementRandomSequencePosition();
        ss.insertItemAtPosition(ssIt, randomPos, false);
        
        if ((i + 1) % 2000 == 0) {
            std::stringstream name;
            name << "After " << (i + 1) << " insertions";
            assertTraversalsMatch(fl, ss, name.str());
            std::cout << "  ✓ " << (i + 1) << " insertions verified (size: " << ss.size() << ")\n";
        }
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    
    assertTraversalsMatch(fl, ss, "Final stress test state");
    
    std::cout << "  Final size: " << ss.size() << "\n";
    std::cout << "  Time taken: " << duration.count() << " ms\n";
    std::cout << "  ✓ STRESS TEST PASSED!\n";
}

void testMultipleSeedsSmallScale() {
    std::cout << "\n=== Test: Multiple seeds with moderate insertions ===\n";
    
    int numSeeds = 10;
    int numInsertions = 1000;
    int initialSize = 50;
    
    for (int i = 0; i < numSeeds; ++i) {
        unsigned int seed = 1000 + i;
        
        size_t requiredCapacity = 1 + initialSize + numInsertions;
        FixedList fl(requiredCapacity);
        fl.initialize(initialSize);
        
        SuperSequence ss(initialSize, 1);
        
        std::mt19937 rng(seed);
        
        for (int j = 0; j < numInsertions; ++j) {
            size_t currentSize = ss.size();
            std::uniform_int_distribution<size_t> dist(0, currentSize - 1);
            size_t insertAfterPos = dist(rng);
            
            auto flIt = fl.begin();
            ++flIt;
            flIt = advanceFixed(flIt, insertAfterPos);
            fl.insertAfter(flIt, false);
            
            auto ssIt = ss.begin();
            ssIt = advanceSuper(ssIt, insertAfterPos);
            ++ssIt;
            size_t randomPos = (j == 0) ? ss.getRandomSequencePosition() : ss.incrementRandomSequencePosition();
            ss.insertItemAtPosition(ssIt, randomPos, false);
        }
        
        std::stringstream name;
        name << "Seed " << seed << " - " << numInsertions << " insertions";
        assertTraversalsMatch(fl, ss, name.str());
        std::cout << "  ✓ Seed " << seed << " passed\n";
    }
    
    std::cout << "  ✓ ALL " << numSeeds << " SEEDS PASSED!\n";
}

void testEdgeCaseInsertions() {
    std::cout << "\n=== Test: Edge case insertion patterns ===\n";
    
    // Test 1: Always insert at beginning
    {
        FixedList fl(1000);
        fl.initialize(10);
        SuperSequence ss(10, 1);
        
        for (int i = 0; i < 100; ++i) {
            auto flIt = fl.begin();
            ++flIt;  // Always at first real position
            fl.insertAfter(flIt, false);
            
            auto ssIt = ss.begin();
            ++ssIt;
            size_t pos = (i == 0) ? ss.getRandomSequencePosition() : ss.incrementRandomSequencePosition();
            ss.insertItemAtPosition(ssIt, pos, false);
        }
        
        assertTraversalsMatch(fl, ss, "Always insert at beginning");
        std::cout << "  ✓ Always insert at beginning\n";
    }
    
    // Test 2: Always insert at end
    {
        FixedList fl(1000);
        fl.initialize(10);
        SuperSequence ss(10, 1);
        
        for (int i = 0; i < 100; ++i) {
            auto flIt = fl.begin();
            ++flIt;
            flIt = advanceFixed(flIt, ss.size() - 1);  // At last position
            fl.insertAfter(flIt, false);
            
            auto ssIt = ss.begin();
            ssIt = advanceSuper(ssIt, ss.size() - 1);
            ++ssIt;
            size_t pos = (i == 0) ? ss.getRandomSequencePosition() : ss.incrementRandomSequencePosition();
            ss.insertItemAtPosition(ssIt, pos, false);
        }
        
        assertTraversalsMatch(fl, ss, "Always insert at end");
        std::cout << "  ✓ Always insert at end\n";
    }
    
    // Test 3: Alternating beginning and end
    {
        FixedList fl(1000);
        fl.initialize(10);
        SuperSequence ss(10, 1);
        
        for (int i = 0; i < 100; ++i) {
            size_t pos = (i % 2 == 0) ? 0 : (ss.size() - 1);
            
            auto flIt = fl.begin();
            ++flIt;
            flIt = advanceFixed(flIt, pos);
            fl.insertAfter(flIt, false);
            
            auto ssIt = ss.begin();
            ssIt = advanceSuper(ssIt, pos);
            ++ssIt;
            size_t randomPos = (i == 0) ? ss.getRandomSequencePosition() : ss.incrementRandomSequencePosition();
            ss.insertItemAtPosition(ssIt, randomPos, false);
        }
        
        assertTraversalsMatch(fl, ss, "Alternating beginning and end");
        std::cout << "  ✓ Alternating beginning and end\n";
    }
}

int main() {
    std::cout << "=========================================\n";
    std::cout << "   Random Insertion Stress Tests        \n";
    std::cout << "=========================================\n";
    
    // Get a time-based seed for truly random tests
    unsigned int timeSeed = static_cast<unsigned int>(
        std::chrono::system_clock::now().time_since_epoch().count()
    );
    
    std::cout << "Time-based seed: " << timeSeed << "\n";
    
    // Progressive difficulty tests
    testRandomInsertions(1000, 100, 12345);
    testRandomInsertions(5000, 50, 54321);
    testRandomInsertions(10000, 100, 99999);
    
    // With references
    testRandomInsertionsWithReferences(5000, 100, 11111);
    
    // Stress test
    testStressTest(10000, 77777);
    
    // Multiple seeds
    testMultipleSeedsSmallScale();
    
    // Edge cases
    testEdgeCaseInsertions();
    
    // Random seed test (will be different each run)
    testRandomInsertions(5000, 50, timeSeed);
    
    std::cout << "\n=========================================\n";
    std::cout << "  All random insertion tests PASSED! ✓  \n";
    std::cout << "=========================================\n";
    std::cout << "\nStatistics:\n";
    std::cout << "- Total insertions tested: 40,000+\n";
    std::cout << "- Different random seeds tested: 10+\n";
    std::cout << "- Edge cases covered: 3\n";
    std::cout << "- With reference tracking: Yes\n";
    std::cout << "\nFixedList robustly matches SuperSequence behavior!\n";
    
    return 0;
}