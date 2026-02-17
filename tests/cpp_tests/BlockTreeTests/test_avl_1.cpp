#include <iostream>
#include <vector>
#include <random>
#include "../../../src/BlockCommon.h"
#include "../../../src/AvlTreeWithRates.h"
#include "../../../src/Event.h"
#include "../../../src/CategorySampler.h"

const size_t MAX_INSERTION_LENGTH = 100;

CategorySampler create_sampler() {
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    return CategorySampler(transMatrix, probs, MAX_INSERTION_LENGTH);
}

void run_test(const std::string& name, bool passed) {
    std::cout << (passed ? "✓ " : "✗ ") << name << "\n";
}

// Test 1: Insert inside insertion part
bool test_insert_in_insertion() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates = {0, 1, 2, 3, 1};
    tree.init_tree(5, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // First insertion at position 2
    Event ev1 = {INSERTION, 2, 3};
    tree.handle_event(ev1, sampler, rng);
    
    // Second insertion inside the first insertion (position 4 is in the insertion)
    Event ev2 = {INSERTION, 4, 2};
    bool success = tree.handle_event(ev2, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 2: Insert at block boundary
bool test_insert_at_boundary() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates = {0, 1, 2, 3, 1};
    tree.init_tree(5, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Create a split at position 2
    Event ev1 = {INSERTION, 2, 1};
    tree.handle_event(ev1, sampler, rng);
    
    // Insert exactly at the boundary (position 3)
    Event ev2 = {INSERTION, 3, 2};
    bool success = tree.handle_event(ev2, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 3: Multiple insertions at same position
bool test_multiple_insertions_same_pos() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates = {0, 1, 2, 3, 1};
    tree.init_tree(5, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Three insertions at position 2
    for (int i = 0; i < 3; i++) {
        Event ev = {INSERTION, 2, 1};
        if (!tree.handle_event(ev, sampler, rng)) return false;
    }
    
    return tree.validate_rate_integrity();
}

// Test 4: Insert at position 0 (start of sequence)
bool test_insert_at_start() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates = {0, 1, 2, 3, 1};
    tree.init_tree(5, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    Event ev = {INSERTION, 0, 3};
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 5: Insert at end then inside that insertion
bool test_insert_at_end_then_inside() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates = {0, 1, 2, 3, 1};
    tree.init_tree(5, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Insert at end
    Event ev1 = {INSERTION, 4, 5};
    tree.handle_event(ev1, sampler, rng);
    
    // Insert inside the insertion we just made
    Event ev2 = {INSERTION, 6, 2};
    bool success = tree.handle_event(ev2, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 6: Large insertion creating multiple blocks
bool test_large_insertion() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(20, 0);
    tree.init_tree(20, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    Event ev = {INSERTION, 10, 15};
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

int main() {
    std::cout << "=== Insertion Edge Case Tests ===\n\n";
    
    int passed = 0;
    int total = 6;
    
    if (test_insert_in_insertion()) passed++;
    run_test("Insert inside insertion", test_insert_in_insertion());
    
    if (test_insert_at_boundary()) passed++;
    run_test("Insert at block boundary", test_insert_at_boundary());
    
    if (test_multiple_insertions_same_pos()) passed++;
    run_test("Multiple insertions at same position", test_multiple_insertions_same_pos());
    
    if (test_insert_at_start()) passed++;
    run_test("Insert at start of sequence", test_insert_at_start());
    
    if (test_insert_at_end_then_inside()) passed++;
    run_test("Insert at end then inside that insertion", test_insert_at_end_then_inside());
    
    if (test_large_insertion()) passed++;
    run_test("Large insertion", test_large_insertion());
    
    std::cout << "\n" << passed << "/" << total << " tests passed\n";
    
    return (passed == total) ? 0 : 1;
}