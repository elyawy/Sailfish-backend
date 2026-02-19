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

// Test 1: Case A - Delete from position 1, only in OP
bool test_delete_case_a() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(20, 0);
    tree.init_tree(20, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    Event ev = {DELETION, 1, 5};  // Delete 5 from position 1 in OP (length=20)
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 2: Case B - Delete entire block (except anchor)
bool test_delete_case_b() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(10, 0);
    tree.init_tree(10, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Add insertion first
    Event ins = {INSERTION, 5, 3};
    tree.handle_event(ins, sampler, rng);
    
    // Delete entire block from position 1 (OP + AP)
    Event ev = {DELETION, 1, 12};
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 3: Case C - Delete from position 1, spanning OP and AP
bool test_delete_case_c() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(20, 0);
    tree.init_tree(20, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Add insertion
    Event ins = {INSERTION, 10, 5};
    tree.handle_event(ins, sampler, rng);
    
    // Delete from position 1, spanning OP and part of AP
    Event ev = {DELETION, 1, 21};  // Delete 20 from OP (positions 1-20), 1 from AP
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 4: Case D - Delete in middle of OP
bool test_delete_case_d() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(20, 0);
    tree.init_tree(20, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    Event ev = {DELETION, 5, 8};  // Delete positions 5-12 in OP
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 5: Case E - Delete ending exactly at OP end
bool test_delete_case_e() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(20, 0);
    tree.init_tree(20, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Add insertion
    Event ins = {INSERTION, 10, 5};
    tree.handle_event(ins, sampler, rng);
    
    // Delete from position 5 to end of OP (position 20)
    Event ev = {DELETION, 5, 15};
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 6: Case F - Delete spanning OP and AP
bool test_delete_case_f() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(20, 0);
    tree.init_tree(20, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Add insertion
    Event ins = {INSERTION, 10, 8};
    tree.handle_event(ins, sampler, rng);
    
    // Delete from middle, spanning OP and AP
    Event ev = {DELETION, 15, 10};  // Delete 5 from OP, 5 from AP
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 7: Multi-block spanning deletion
bool test_multi_block_deletion() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(30, 0);
    tree.init_tree(30, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Create multiple blocks by inserting at different positions
    Event ins1 = {INSERTION, 10, 3};
    tree.handle_event(ins1, sampler, rng);
    Event ins2 = {INSERTION, 20, 3};
    tree.handle_event(ins2, sampler, rng);
    
    // Delete spanning multiple blocks
    Event ev = {DELETION, 8, 20};
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 8: Delete from insertion only
bool test_delete_from_insertion() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(20, 0);
    tree.init_tree(20, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Add large insertion
    Event ins = {INSERTION, 10, 15};
    tree.handle_event(ins, sampler, rng);
    
    // Delete only from the inserted part
    Event ev = {DELETION, 22, 5};  // Position 22 is 2 into the AP
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 9: Delete from first block starting at position 1
bool test_delete_first_block() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(10, 0);
    tree.init_tree(10, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    Event ev = {DELETION, 1, 3};
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

// Test 10: Sequential deletions
bool test_sequential_deletions() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(30, 0);
    tree.init_tree(30, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Multiple deletions at same position (position 5)
    for (int i = 0; i < 3; i++) {
        Event ev = {DELETION, 5, 2};
        if (!tree.handle_event(ev, sampler, rng)) return false;
    }
    
    return tree.validate_rate_integrity();
}

// Test 11: Delete entire middle block
bool test_delete_middle_block() {
    avl_array_with_rates<size_t, size_t, 100> tree;
    std::vector<size_t> parentRates(30, 0);
    tree.init_tree(30, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(42);
    
    // Create 3 blocks
    Event ins1 = {INSERTION, 10, 2};
    tree.handle_event(ins1, sampler, rng);
    Event ins2 = {INSERTION, 20, 2};
    tree.handle_event(ins2, sampler, rng);
    
    // Delete entire middle block and more
    Event ev = {DELETION, 8, 18};
    bool success = tree.handle_event(ev, sampler, rng);
    
    return success && tree.validate_rate_integrity();
}

int main() {
    std::cout << "=== Deletion Edge Case Tests ===\n\n";
    
    int passed = 0;
    int total = 11;
    
    if (test_delete_case_a()) passed++;
    run_test("Case A: Delete from position 1 (OP only)", test_delete_case_a());
    
    if (test_delete_case_b()) passed++;
    run_test("Case B: Delete entire block from position 1", test_delete_case_b());
    
    if (test_delete_case_c()) passed++;
    run_test("Case C: Delete from position 1 (spanning OP+AP)", test_delete_case_c());
    
    if (test_delete_case_d()) passed++;
    run_test("Case D: Delete in middle of OP", test_delete_case_d());
    
    if (test_delete_case_e()) passed++;
    run_test("Case E: Delete ending at OP end", test_delete_case_e());
    
    if (test_delete_case_f()) passed++;
    run_test("Case F: Delete spanning OP and AP", test_delete_case_f());
    
    if (test_multi_block_deletion()) passed++;
    run_test("Multi-block spanning deletion", test_multi_block_deletion());
    
    if (test_delete_from_insertion()) passed++;
    run_test("Delete from insertion only", test_delete_from_insertion());
    
    if (test_delete_first_block()) passed++;
    run_test("Delete from first block at position 1", test_delete_first_block());
    
    if (test_sequential_deletions()) passed++;
    run_test("Sequential deletions at same position", test_sequential_deletions());
    
    if (test_delete_middle_block()) passed++;
    run_test("Delete entire middle block", test_delete_middle_block());
    
    std::cout << "\n" << passed << "/" << total << " tests passed\n";
    
    return (passed == total) ? 0 : 1;
}