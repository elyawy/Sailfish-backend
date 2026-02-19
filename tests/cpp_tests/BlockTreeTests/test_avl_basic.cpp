#include <iostream>
#include <vector>
#include <random>
#include "../../../src/BlockCommon.h"
#include "../../../src/AvlTreeWithRates.h"
#include "../../../src/Event.h"
#include "../../../src/CategorySampler.h"

const size_t MAX_INSERTION_LENGTH = 100;

// Helper to create a simple sampler
CategorySampler create_simple_sampler() {
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    return CategorySampler(transMatrix, probs, MAX_INSERTION_LENGTH);
}

// Test 1: Insert 1 element into a block of length 1
bool test_insert_into_length_1() {
    std::cout << "\n========================================\n";
    std::cout << "TEST 1: Insert 1 element into block of length 1\n";
    std::cout << "Expected: Block grows from [length=1, insertion=0] to [length=1, insertion=1]\n";
    std::cout << "Expected: rateCategories grows from size 0 to size 1\n";
    std::cout << "========================================\n";
    
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates = {0}; // Single parent rate
    tree.init_tree(1, parentRates);
    
    std::cout << "Initial tree:\n" << tree.print_avl() << "\n";
    
    // Verify initial state
    auto it = tree.begin();
    std::cout << "Initial block: length=" << (*it).length 
              << " insertion=" << (*it).insertion
              << " rateCategories.size()=" << (*it).rateCategories.size() << "\n";
    
    if ((*it).length != 1 || (*it).insertion != 0 || (*it).rateCategories.size() != 0) {
        std::cout << "FAIL: Initial state incorrect\n";
        return false;
    }
    
    // Insert 1 element at position 0 (after position 0 with +1 semantics)
    CategorySampler sampler = create_simple_sampler();
    std::mt19937_64 rng(42);
    Event ev = {INSERTION, 0, 1};
    
    std::cout << "\nInserting " << ev.length << " element(s) at position " << ev.position << "\n";
    bool success = tree.handle_event(ev, sampler, rng);
    
    std::cout << "Tree after insertion:\n" << tree.print_avl() << "\n";
    
    // Verify final state
    it = tree.begin();
    std::cout << "Final block: length=" << (*it).length 
              << " insertion=" << (*it).insertion
              << " rateCategories.size()=" << (*it).rateCategories.size() << "\n";
    
    if (!success) {
        std::cout << "FAIL: Insertion returned false\n";
        return false;
    }
    
    if ((*it).length != 1 || (*it).insertion != 1 || (*it).rateCategories.size() != 1) {
        std::cout << "FAIL: Final state incorrect\n";
        std::cout << "  Expected: length=1, insertion=1, rateCategories.size()=1\n";
        std::cout << "  Got: length=" << (*it).length 
                  << ", insertion=" << (*it).insertion
                  << ", rateCategories.size()=" << (*it).rateCategories.size() << "\n";
        return false;
    }
    
    std::cout << "PASS\n";
    return true;
}

// Test 2: Insert 1 element in the middle of a block of length 5
bool test_insert_in_middle_of_5() {
    std::cout << "\n========================================\n";
    std::cout << "TEST 2: Insert 1 element in middle of block of length 5\n";
    std::cout << "Expected: Block splits into two blocks\n";
    std::cout << "========================================\n";
    
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates = {0, 1, 2, 3, 1}; // 5 parent rates
    tree.init_tree(5, parentRates);
    
    std::cout << "Initial tree:\n" << tree.print_avl() << "\n";
    
    // Insert at position 2 (middle of block)
    CategorySampler sampler = create_simple_sampler();
    std::mt19937_64 rng(42);
    Event ev = {INSERTION, 2, 1};
    
    std::cout << "\nInserting " << ev.length << " element(s) at position " << ev.position << "\n";
    std::cout << "Expected split: [0|3|0] and [2+insertion(1)|2|0]\n";
    
    bool success = tree.handle_event(ev, sampler, rng);
    
    std::cout << "Tree after insertion:\n" << tree.print_avl() << "\n";
    
    if (!success) {
        std::cout << "FAIL: Insertion returned false\n";
        return false;
    }
    
    // Count blocks
    int block_count = 0;
    for (auto it = tree.begin(); it != tree.end(); ++it) {
        block_count++;
        std::cout << "Block " << block_count << ": key=" << it.key()
                  << " length=" << (*it).length
                  << " insertion=" << (*it).insertion
                  << " rateCategories.size()=" << (*it).rateCategories.size() << "\n";
    }
    
    if (block_count != 2) {
        std::cout << "FAIL: Expected 2 blocks, got " << block_count << "\n";
        return false;
    }
    
    bool valid = tree.validate_rate_integrity();
    if (!valid) {
        std::cout << "FAIL: Rate integrity check failed\n";
        return false;
    }
    
    std::cout << "PASS\n";
    return true;
}

// Test 3: Insert 1 element at the end of a block of length 5
bool test_insert_at_end_of_5() {
    std::cout << "\n========================================\n";
    std::cout << "TEST 3: Insert 1 element at end of block of length 5\n";
    std::cout << "Expected: Single block with insertion added\n";
    std::cout << "========================================\n";
    
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates = {0, 1, 2, 3, 1}; // 5 parent rates
    tree.init_tree(5, parentRates);
    
    std::cout << "Initial tree:\n" << tree.print_avl() << "\n";
    
    // Insert at position 4 (last position, index 4 = after position 4 with +1 semantics)
    CategorySampler sampler = create_simple_sampler();
    std::mt19937_64 rng(42);
    Event ev = {INSERTION, 4, 1};
    
    std::cout << "\nInserting " << ev.length << " element(s) at position " << ev.position << "\n";
    std::cout << "Expected: [0|5|1] with rateCategories.size()=1\n";
    
    bool success = tree.handle_event(ev, sampler, rng);
    
    std::cout << "Tree after insertion:\n" << tree.print_avl() << "\n";
    
    if (!success) {
        std::cout << "FAIL: Insertion returned false\n";
        return false;
    }
    
    // Verify single block
    int block_count = 0;
    auto it = tree.begin();
    for (; it != tree.end(); ++it) {
        block_count++;
        std::cout << "Block " << block_count << ": key=" << it.key()
                  << " length=" << (*it).length
                  << " insertion=" << (*it).insertion
                  << " rateCategories.size()=" << (*it).rateCategories.size() << "\n";
    }
    
    if (block_count != 1) {
        std::cout << "FAIL: Expected 1 block, got " << block_count << "\n";
        return false;
    }
    
    it = tree.begin();
    if ((*it).length != 5 || (*it).insertion != 1 || (*it).rateCategories.size() != 1) {
        std::cout << "FAIL: Final state incorrect\n";
        std::cout << "  Expected: length=5, insertion=1, rateCategories.size()=1\n";
        std::cout << "  Got: length=" << (*it).length 
                  << ", insertion=" << (*it).insertion
                  << ", rateCategories.size()=" << (*it).rateCategories.size() << "\n";
        return false;
    }
    
    bool valid = tree.validate_rate_integrity();
    if (!valid) {
        std::cout << "FAIL: Rate integrity check failed\n";
        return false;
    }
    
    std::cout << "PASS\n";
    return true;
}

int main() {
    std::cout << "========================================\n";
    std::cout << "Basic AvlTreeWithRates Tests\n";
    std::cout << "========================================\n";
    
    int passed = 0;
    int total = 3;
    
    if (test_insert_into_length_1()) passed++;
    if (test_insert_in_middle_of_5()) passed++;
    if (test_insert_at_end_of_5()) passed++;
    
    std::cout << "\n========================================\n";
    std::cout << "RESULTS: " << passed << "/" << total << " tests passed\n";
    std::cout << "========================================\n";
    
    return (passed == total) ? 0 : 1;
}