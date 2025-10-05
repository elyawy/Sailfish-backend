#include <iostream>
#include <vector>
#include <cassert>
#include <sstream>
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
    } else {
        std::cout << "PASSED: " << testName << "\n";
        std::cout << "  Matching traversal: " << fixedTraversal << "\n";
    }
}

void testInitialization() {
    std::cout << "\n=== Test: Initialization ===\n";
    
    FixedList fl(20);
    fl.initialize(10);
    
    SuperSequence ss(10, 1);  // sequenceSize=10, numSequences=1
    
    assert(fl.size() == 11);  // 1 initial (index 0) + 10
    assert(ss.size() == 10);  // 10 elements with positions 1-10
    
    assertTraversalsMatch(fl, ss, "After initialization with 10 elements");
}

void testInsertAfterFirst() {
    std::cout << "\n=== Test: Insert after first element ===\n";
    
    FixedList fl(20);
    fl.initialize(5);  // Creates: 0-1-2-3-4-5 (skipping 0: 1-2-3-4-5)
    
    SuperSequence ss(5, 1);  // Creates: 1-2-3-4-5
    
    // Insert after first element in both
    auto flIt = fl.begin();
    ++flIt;  // Move to element at index 1 (first real element after phantom)
    fl.insertAfter(flIt, false);
    
    auto ssIt = ss.begin();  // Points to position 1
    ++ssIt;  // Move past position 1
    size_t randomPos = ss.getRandomSequencePosition();
    ss.insertItemAtPosition(ssIt, randomPos, false);
    
    assertTraversalsMatch(fl, ss, "After inserting after first element");
}

void testInsertAfterMiddle() {
    std::cout << "\n=== Test: Insert after middle element ===\n";
    
    FixedList fl(20);
    fl.initialize(5);
    
    SuperSequence ss(5, 1);
    
    // Insert after position 3 in both
    auto flIt = fl.begin();
    ++flIt;  // Skip phantom
    flIt = advanceFixed(flIt, 2);  // Now at index 3 (position 3)
    fl.insertAfter(flIt, false);
    
    auto ssIt = ss.begin();
    ssIt = advanceSuper(ssIt, 2);  // At position 3
    ++ssIt;  // Move past position 3
    size_t randomPos = ss.getRandomSequencePosition();
    ss.insertItemAtPosition(ssIt, randomPos, false);
    
    assertTraversalsMatch(fl, ss, "After inserting after middle element");
}

void testInsertAfterLast() {
    std::cout << "\n=== Test: Insert after last element ===\n";
    
    FixedList fl(20);
    fl.initialize(5);
    
    SuperSequence ss(5, 1);
    
    // Insert after last element in both
    auto flIt = fl.begin();
    ++flIt;  // Skip phantom
    flIt = advanceFixed(flIt, 4);  // Now at index 5 (last element)
    fl.insertAfter(flIt, false);
    
    auto ssIt = ss.begin();
    ssIt = advanceSuper(ssIt, 4);  // At position 5 (last)
    ++ssIt;  // Move to end()
    size_t randomPos = ss.getRandomSequencePosition();
    ss.insertItemAtPosition(ssIt, randomPos, false);
    
    assertTraversalsMatch(fl, ss, "After inserting after last element");
}

void testMultipleInsertions() {
    std::cout << "\n=== Test: Multiple insertions at various positions ===\n";
    
    FixedList fl(20);
    fl.initialize(5);
    
    SuperSequence ss(5, 1);
    
    // Insert #1: After position 1
    auto flIt = fl.begin();
    ++flIt;  // Skip phantom, now at index 1
    fl.insertAfter(flIt, false);
    
    auto ssIt = ss.begin();  // At position 1
    ++ssIt;  // Move past position 1
    size_t pos1 = ss.getRandomSequencePosition();
    ss.insertItemAtPosition(ssIt, pos1, false);
    
    assertTraversalsMatch(fl, ss, "After 1st insertion");
    
    // Insert #2: After position 3 (in the new sequence)
    flIt = fl.begin();
    ++flIt;  // Skip phantom
    flIt = advanceFixed(flIt, 2);  // At position 3
    fl.insertAfter(flIt, false);
    
    ssIt = ss.begin();
    ssIt = advanceSuper(ssIt, 2);  // At position 3
    ++ssIt;
    size_t pos2 = ss.incrementRandomSequencePosition();
    ss.insertItemAtPosition(ssIt, pos2, false);
    
    assertTraversalsMatch(fl, ss, "After 2nd insertion");
    
    // Insert #3: After last element
    flIt = fl.begin();
    ++flIt;  // Skip phantom
    flIt = advanceFixed(flIt, fl.size() - 2);  // At last real element
    fl.insertAfter(flIt, false);
    
    ssIt = ss.begin();
    ssIt = advanceSuper(ssIt, ss.size() - 1);  // At last element
    ++ssIt;
    size_t pos3 = ss.incrementRandomSequencePosition();
    ss.insertItemAtPosition(ssIt, pos3, false);
    
    assertTraversalsMatch(fl, ss, "After 3rd insertion");
}

void testSequentialInsertions() {
    std::cout << "\n=== Test: Sequential insertions building up chain ===\n";
    
    FixedList fl(20);
    fl.initialize(3);
    
    SuperSequence ss(3, 1);
    
    assertTraversalsMatch(fl, ss, "Initial state");
    
    // Build up: Insert after first element three times
    for (int i = 0; i < 3; i++) {
        auto flIt = fl.begin();
        ++flIt;  // Skip phantom
        flIt = advanceFixed(flIt, i);  // Move to position where we want to insert
        fl.insertAfter(flIt, false);
        
        auto ssIt = ss.begin();
        ssIt = advanceSuper(ssIt, i);
        ++ssIt;
        size_t pos = (i == 0) ? ss.getRandomSequencePosition() : ss.incrementRandomSequencePosition();
        ss.insertItemAtPosition(ssIt, pos, false);
        
        std::stringstream name;
        name << "After insertion #" << (i + 1);
        assertTraversalsMatch(fl, ss, name.str());
    }
}

void testReferenceTracking() {
    std::cout << "\n=== Test: Reference tracking (isColumn) ===\n";
    
    FixedList fl(20);
    fl.initialize(5);
    
    SuperSequence ss(5, 1);
    
    // Reference positions 1, 3, 5 (0-indexed: positions 0, 2, 4 after phantom)
    auto flIt = fl.begin();
    ++flIt;  // Skip phantom, now at index 1
    size_t idx1 = *flIt;
    fl.referencePosition(flIt);
    
    flIt = fl.begin();
    ++flIt;  // Skip phantom
    flIt = advanceFixed(flIt, 2);  // At index 3
    size_t idx3 = *flIt;
    fl.referencePosition(flIt);
    
    flIt = fl.begin();
    ++flIt;  // Skip phantom
    flIt = advanceFixed(flIt, 4);  // At index 5
    size_t idx5 = *flIt;
    fl.referencePosition(flIt);
    
    // Same for SuperSequence
    auto ssIt = ss.begin();
    ss.referencePosition(ssIt);
    
    ssIt = advanceSuper(ss.begin(), 2);
    ss.referencePosition(ssIt);
    
    ssIt = advanceSuper(ss.begin(), 4);
    ss.referencePosition(ssIt);
    
    // Check isColumn flags match
    assert(fl.getIsColumn(idx1) == true);
    assert(fl.getIsColumn(idx3) == true);
    assert(fl.getIsColumn(idx5) == true);
    
    ssIt = ss.begin();
    assert(ssIt->isColumn == true);
    ssIt = advanceSuper(ss.begin(), 2);
    assert(ssIt->isColumn == true);
    ssIt = advanceSuper(ss.begin(), 4);
    assert(ssIt->isColumn == true);
    
    assertTraversalsMatch(fl, ss, "After setting reference flags");
    
    std::cout << "  Reference tracking matches correctly\n";
}

void testAbsolutePositions() {
    std::cout << "\n=== Test: Absolute positions after setAbsolutePositions ===\n";
    
    FixedList fl(20);
    fl.initialize(5);
    
    SuperSequence ss(5, 1);
    
    // Reference some positions (1 and 4)
    auto flIt = fl.begin();
    ++flIt;  // Skip phantom, at index 1
    size_t idx1 = *flIt;
    fl.referencePosition(flIt);
    
    flIt = fl.begin();
    ++flIt;  // Skip phantom
    flIt = advanceFixed(flIt, 3);  // At index 4
    size_t idx4 = *flIt;
    fl.referencePosition(flIt);
    
    auto ssIt = ss.begin();
    ss.referencePosition(ssIt);
    
    ssIt = advanceSuper(ss.begin(), 3);
    ss.referencePosition(ssIt);
    
    // Set absolute positions
    fl.setAbsolutePositions();
    ss.setAbsolutePositions();
    
    assert(fl.getMsaSequenceLength() == 2);
    assert(ss.getMsaSequenceLength() == 2);
    std::cout << "  MSA lengths match: " << fl.getMsaSequenceLength() << "\n";
    
    // Check absolute positions match
    assert(fl.getAbsolutePosition(idx1) == 0);
    assert(fl.getAbsolutePosition(idx4) == 1);
    
    ssIt = ss.begin();
    assert(ssIt->absolutePosition == 0);
    ssIt = advanceSuper(ss.begin(), 3);
    assert(ssIt->absolutePosition == 1);
    
    std::cout << "  Absolute positions match correctly\n";
    
    assertTraversalsMatch(fl, ss, "After setAbsolutePositions");
}

void testInsertWithReference() {
    std::cout << "\n=== Test: Insert with reference flag ===\n";
    
    FixedList fl(20);
    fl.initialize(3);
    
    SuperSequence ss(3, 1);
    
    // Insert after first with reference=true
    auto flIt = fl.begin();
    ++flIt;  // Skip phantom, at index 1
    fl.insertAfter(flIt, true);
    
    auto ssIt = ss.begin();
    ++ssIt;
    size_t randomPos = ss.getRandomSequencePosition();
    ss.insertItemAtPosition(ssIt, randomPos, true);
    
    assertTraversalsMatch(fl, ss, "After insert with reference=true");
    
    // Set absolute positions and check MSA length
    fl.setAbsolutePositions();
    ss.setAbsolutePositions();
    
    assert(fl.getMsaSequenceLength() == 1);
    assert(ss.getMsaSequenceLength() == 1);
    
    std::cout << "  MSA lengths match: " << fl.getMsaSequenceLength() << "\n";
}

void testComplexSequence() {
    std::cout << "\n=== Test: Complex sequence of operations ===\n";
    
    FixedList fl(30);
    fl.initialize(5);
    
    SuperSequence ss(5, 2);  // numSequences=2 for testing
    
    assertTraversalsMatch(fl, ss, "Initial state");
    
    // Operation 1: Insert after position 2
    auto flIt = fl.begin();
    ++flIt;  // Skip phantom
    ++flIt;  // At index 2
    fl.insertAfter(flIt, false);
    
    auto ssIt = ss.begin();
    ++ssIt;  // At position 2
    ++ssIt;  // Move past position 2
    size_t pos = ss.getRandomSequencePosition();
    ss.insertItemAtPosition(ssIt, pos, false);
    
    assertTraversalsMatch(fl, ss, "After insert at position 2");
    
    // Operation 2: Reference position 3
    flIt = fl.begin();
    ++flIt;  // Skip phantom
    flIt = advanceFixed(flIt, 2);  // At index 3
    fl.referencePosition(flIt);
    
    ssIt = ss.begin();
    ssIt = advanceSuper(ssIt, 2);  // At position 3
    ss.referencePosition(ssIt);
    
    assertTraversalsMatch(fl, ss, "After referencing position 3");
    
    // Operation 3: Insert after position 4 with reference
    flIt = fl.begin();
    ++flIt;  // Skip phantom
    flIt = advanceFixed(flIt, 3);  // At index 4
    fl.insertAfter(flIt, true);
    
    ssIt = ss.begin();
    ssIt = advanceSuper(ssIt, 3);  // At position 4
    ++ssIt;
    pos = ss.incrementRandomSequencePosition();
    ss.insertItemAtPosition(ssIt, pos, true);
    
    assertTraversalsMatch(fl, ss, "After insert with reference at position 4");
    
    // Operation 4: Set absolute positions
    fl.setAbsolutePositions();
    ss.setAbsolutePositions();
    
    assert(fl.getMsaSequenceLength() == ss.getMsaSequenceLength());
    std::cout << "  Final MSA length: " << fl.getMsaSequenceLength() << "\n";
    
    assertTraversalsMatch(fl, ss, "Final state after setAbsolutePositions");
}

void testEdgeCases() {
    std::cout << "\n=== Test: Edge cases ===\n";
    
    // Test with size 1
    FixedList fl1(10);
    fl1.initialize(1);
    SuperSequence ss1(1, 1);
    assertTraversalsMatch(fl1, ss1, "Size 1 initialization");
    
    // Test with size 2
    FixedList fl2(10);
    fl2.initialize(2);
    SuperSequence ss2(2, 1);
    assertTraversalsMatch(fl2, ss2, "Size 2 initialization");
    
    // Insert after only element
    auto flIt = fl1.begin();
    ++flIt;
    fl1.insertAfter(flIt, false);
    
    auto ssIt = ss1.begin();
    ++ssIt;
    ss1.insertItemAtPosition(ssIt, ss1.getRandomSequencePosition(), false);
    
    assertTraversalsMatch(fl1, ss1, "Insert after single element");
}

int main() {
    std::cout << "=========================================\n";
    std::cout << " FixedList vs SuperSequence Comparison  \n";
    std::cout << "=========================================\n";
    
    testInitialization();
    testInsertAfterFirst();
    testInsertAfterMiddle();
    testInsertAfterLast();
    testMultipleInsertions();
    testSequentialInsertions();
    testReferenceTracking();
    testAbsolutePositions();
    testInsertWithReference();
    testComplexSequence();
    testEdgeCases();
    
    std::cout << "\n=========================================\n";
    std::cout << "  All comparison tests PASSED! ✓        \n";
    std::cout << "=========================================\n";
    std::cout << "\nKey findings:\n";
    std::cout << "1. Traversal sequences match perfectly (1-based positions)\n";
    std::cout << "2. Insertion behavior is equivalent\n";
    std::cout << "3. Reference tracking (isColumn) matches\n";
    std::cout << "4. MSA length computation is identical\n";
    std::cout << "5. Absolute position assignment matches\n";
    std::cout << "\nFixedList successfully replicates SuperSequence behavior!\n";
    
    return 0;
}