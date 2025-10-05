#include <iostream>
#include <vector>
#include <cassert>
#include <sstream>
#include "../../src/FixedList.h"

// Helper function to advance iterator manually
FixedList::iterator advance(FixedList::iterator it, int n) {
    for (int i = 0; i < n; ++i) {
        ++it;
    }
    return it;
}

// Helper function to get the traversal as a string
std::string getTraversal(FixedList& fl) {
    std::stringstream ss;
    bool first = true;
    for (auto it = fl.begin(); it != fl.end(); ++it) {
        if (!first) ss << "-";
        ss << *it;
        first = false;
    }
    return ss.str();
}

// Helper function to verify traversal matches expected
void assertTraversal(FixedList& fl, const std::string& expected, const std::string& test_name) {
    std::string actual = getTraversal(fl);
    if (actual != expected) {
        std::cout << "FAILED: " << test_name << "\n";
        std::cout << "  Expected: " << expected << "\n";
        std::cout << "  Actual:   " << actual << "\n";
        assert(false);
    } else {
        std::cout << "PASSED: " << test_name << "\n";
        std::cout << "  Traversal: " << actual << "\n";
    }
}

void testInitialization() {
    std::cout << "\n=== Test: Initialization ===\n";
    
    FixedList fl(20);
    fl.initialize(10);
    
    assert(fl.size() == 11); // 1 initial + 10 inserted
    assertTraversal(fl, "0-1-2-3-4-5-6-7-8-9-10", "Initial 10 elements");
    
    // Check MSA length starts at 0 (no references yet)
    assert(fl.getMsaSequenceLength() == SIZE_MAX); // INVALID
}

void testInsertAfterPosition0() {
    std::cout << "\n=== Test: Insert after position 0 (after first element) ===\n";
    
    FixedList fl(20);
    fl.initialize(10);
    
    auto it = fl.begin();  // Position 0 (first element)
    fl.insertAfter(it, false);
    
    assert(fl.size() == 12); // 1 initial + 10 + 1 inserted
    assertTraversal(fl, "0-11-1-2-3-4-5-6-7-8-9-10", "Insert after position 0");
}

void testInsertAfterPosition1() {
    std::cout << "\n=== Test: Insert after position 1 ===\n";
    
    FixedList fl(20);
    fl.initialize(10);
    
    auto it = fl.begin();
    ++it;  // Position 1
    fl.insertAfter(it, false);
    
    assert(fl.size() == 12);
    assertTraversal(fl, "0-1-11-2-3-4-5-6-7-8-9-10", "Insert after position 1");
}

void testInsertAfterPosition5() {
    std::cout << "\n=== Test: Insert after position 5 (middle) ===\n";
    
    FixedList fl(20);
    fl.initialize(10);
    
    auto it = advance(fl.begin(), 5);  // Position 5
    fl.insertAfter(it, false);
    
    assert(fl.size() == 12);
    assertTraversal(fl, "0-1-2-3-4-5-11-6-7-8-9-10", "Insert after position 5");
}

void testInsertAfterPosition9() {
    std::cout << "\n=== Test: Insert after position 9 (last element) ===\n";
    
    FixedList fl(20);
    fl.initialize(10);
    
    auto it = advance(fl.begin(), 9);  // Position 9
    fl.insertAfter(it, false);
    
    assert(fl.size() == 12);
    assertTraversal(fl, "0-1-2-3-4-5-6-7-8-9-11-10", "Insert after position 9");
}

void testAllInsertionPositions() {
    std::cout << "\n=== Test: All 11 insertion positions (after each element) ===\n";
    
    std::vector<std::string> expected = {
        "0-11-1-2-3-4-5-6-7-8-9-10",     // after pos 0
        "0-1-11-2-3-4-5-6-7-8-9-10",     // after pos 1
        "0-1-2-11-3-4-5-6-7-8-9-10",     // after pos 2
        "0-1-2-3-11-4-5-6-7-8-9-10",     // after pos 3
        "0-1-2-3-4-11-5-6-7-8-9-10",     // after pos 4
        "0-1-2-3-4-5-11-6-7-8-9-10",     // after pos 5
        "0-1-2-3-4-5-6-11-7-8-9-10",     // after pos 6
        "0-1-2-3-4-5-6-7-11-8-9-10",     // after pos 7
        "0-1-2-3-4-5-6-7-8-11-9-10",     // after pos 8
        "0-1-2-3-4-5-6-7-8-9-11-10",     // after pos 9
        "0-1-2-3-4-5-6-7-8-9-10-11"      // after pos 10
    };
    
    for (int insert_after_pos = 0; insert_after_pos < 11; insert_after_pos++) {
        FixedList fl(20);
        fl.initialize(10);
        
        auto it = advance(fl.begin(), insert_after_pos);
        
        fl.insertAfter(it, false);
        
        std::stringstream test_name;
        test_name << "Insert after position " << insert_after_pos;
        assertTraversal(fl, expected[insert_after_pos], test_name.str());
    }
}

void testMultipleInsertions() {
    std::cout << "\n=== Test: Multiple insertions ===\n";
    
    FixedList fl(20);
    fl.initialize(10);
    
    // Insert after position 5
    auto it = advance(fl.begin(), 5);
    fl.insertAfter(it, false);
    assertTraversal(fl, "0-1-2-3-4-5-11-6-7-8-9-10", "After first insertion");
    
    // Insert after position 2
    it = advance(fl.begin(), 2);
    fl.insertAfter(it, false);
    assertTraversal(fl, "0-1-2-12-3-4-5-11-6-7-8-9-10", "After second insertion");
    
    // Insert after last element
    it = advance(fl.begin(), fl.size() - 1);
    fl.insertAfter(it, false);
    assertTraversal(fl, "0-1-2-12-3-4-5-11-6-7-8-9-10-13", "After third insertion");
    
    assert(fl.size() == 14); // 1 initial + 10 + 3 inserted
}

void testReferenceTracking() {
    std::cout << "\n=== Test: Reference tracking (MSA length) ===\n";
    
    FixedList fl(20);
    fl.initialize(10);
    
    assert(fl.getMsaSequenceLength() == SIZE_MAX); // INVALID initially
    std::cout << "Initial MSA length: INVALID\n";
    
    // Reference an existing position
    auto it = fl.begin();
    fl.referencePosition(it);
    std::cout << "After first reference (should not update until setAbsolutePositions)\n";
    
    // Reference another position
    it = advance(fl.begin(), 3);
    fl.referencePosition(it);
    
    // Set absolute positions to compute MSA length
    fl.setAbsolutePositions();
    
    assert(fl.getMsaSequenceLength() == 2);
    std::cout << "After setAbsolutePositions: 2\n";
}

void testCapacity() {
    std::cout << "\n=== Test: Filling to capacity ===\n";
    
    FixedList fl(16);  // Exact capacity for 1 initial + 10 initialize + 5 insertions
    fl.initialize(10);
    
    // Insert 5 elements
    for (int i = 0; i < 5; i++) {
        auto it = advance(fl.begin(), i * 2);  // Insert after positions 0, 2, 4, 6, 8
        fl.insertAfter(it, false);
    }
    
    assert(fl.size() == 16);
    assertTraversal(fl, "0-11-1-12-2-13-3-14-4-15-5-6-7-8-9-10", "Filled to capacity");
    
    std::cout << "Successfully filled to exact capacity\n";
}

void testAbsolutePositions() {
    std::cout << "\n=== Test: Absolute positions after setAbsolutePositions ===\n";
    
    FixedList fl(20);
    fl.initialize(5);
    
    // Reference some positions (indices 0, 2, 4)
    auto it = fl.begin();
    size_t idx0 = *it;
    fl.referencePosition(it);
    
    it = advance(fl.begin(), 2);
    size_t idx2 = *it;
    fl.referencePosition(it);
    
    it = advance(fl.begin(), 4);
    size_t idx4 = *it;
    fl.referencePosition(it);
    
    // Set absolute positions
    fl.setAbsolutePositions();
    
    assert(fl.getMsaSequenceLength() == 3);
    
    // Check absolute positions
    assert(fl.getIsColumn(idx0));
    assert(fl.getAbsolutePosition(idx0) == 0);
    std::cout << "Index " << idx0 << " -> absolute 0\n";
    
    assert(fl.getIsColumn(idx2));
    assert(fl.getAbsolutePosition(idx2) == 1);
    std::cout << "Index " << idx2 << " -> absolute 1\n";
    
    assert(fl.getIsColumn(idx4));
    assert(fl.getAbsolutePosition(idx4) == 2);
    std::cout << "Index " << idx4 << " -> absolute 2\n";
}

void testSequentialInsertions() {
    std::cout << "\n=== Test: Sequential insertions (building up) ===\n";
    
    FixedList fl(20);
    fl.initialize(3);
    
    assertTraversal(fl, "0-1-2-3", "Initial 3 elements");
    
    // Insert after first element repeatedly
    auto it = fl.begin();
    fl.insertAfter(it, false);
    assertTraversal(fl, "0-4-1-2-3", "After 1st insertion");
    
    it = advance(fl.begin(), 1);  // Now at the newly inserted element
    fl.insertAfter(it, false);
    assertTraversal(fl, "0-4-5-1-2-3", "After 2nd insertion");
    
    it = advance(fl.begin(), 2);
    fl.insertAfter(it, false);
    assertTraversal(fl, "0-4-5-6-1-2-3", "After 3rd insertion");
}

int main() {
    std::cout << "====================================\n";
    std::cout << "    FixedList Comprehensive Tests   \n";
    std::cout << "====================================\n";
    
    testInitialization();
    testInsertAfterPosition0();
    testInsertAfterPosition1();
    testInsertAfterPosition5();
    testInsertAfterPosition9();
    testAllInsertionPositions();
    testMultipleInsertions();
    testReferenceTracking();
    testCapacity();
    testAbsolutePositions();
    testSequentialInsertions();
    
    std::cout << "\n====================================\n";
    std::cout << "   All FixedList tests PASSED! ✓   \n";
    std::cout << "====================================\n";
    
    return 0;
}