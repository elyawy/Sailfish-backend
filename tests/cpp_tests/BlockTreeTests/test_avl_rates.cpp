#include <iostream>
#include <vector>
#include <random>
#include "../../../src/BlockCommon.h"
#include "../../../src/AvlTreeWithRates.h"
#include "../../../src/Event.h"
#include "../../../src/CategorySampler.h"

const size_t MAX_INSERTION_LENGTH = 100;

// Helper function to print test header
void print_test_header(const std::string& test_name) {
    std::cout << "\n========================================\n";
    std::cout << "TEST: " << test_name << "\n";
    std::cout << "========================================\n";
}

// Helper function to print event
void print_event(const Event& ev) {
    std::cout << "\nEvent: " << (ev.type == INSERTION ? "INSERTION" : "DELETION")
              << " at position " << ev.position << ", size " << ev.length << "\n";
}

// Test 1: Single insertion in original part
bool test_single_insertion_in_original() {
    print_test_header("Single Insertion in Original Part");
    
    avl_array_with_rates<size_t, size_t, 1000> tree;
    std::vector<size_t> parentRates(100, 0);
    
    tree.init_tree(100, parentRates);
    
    std::cout << "Initial tree:\n" << tree.print_avl() << "\n";
    
    // Create simple category sampler with 4 equal categories
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    CategorySampler sampler(transMatrix, probs, MAX_INSERTION_LENGTH);
    
    std::mt19937_64 rng(42);
    
    // Insert 5 positions at position 50
    Event ev = {INSERTION, 50, 5};
    print_event(ev);
    
    bool success = tree.handle_event(ev, sampler, rng);
    std::cout << "Event success: " << (success ? "YES" : "NO") << "\n";
    std::cout << "Tree after insertion:\n" << tree.print_avl() << "\n";
    
    bool valid = tree.validate_rate_integrity();
    std::cout << "Rate integrity: " << (valid ? "PASS" : "FAIL") << "\n";
    
    return success && valid;
}

// Test 2: Single insertion in added part
bool test_single_insertion_in_added() {
    print_test_header("Single Insertion in Added Part");
    
    avl_array_with_rates<size_t, size_t, 1000> tree;
    std::vector<size_t> parentRates(100, 0);
    tree.init_tree(100, parentRates);
    
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    CategorySampler sampler(transMatrix, probs, MAX_INSERTION_LENGTH);
    std::mt19937_64 rng(42);
    
    // First insertion at position 50 to create a split
    Event ev1 = {INSERTION, 50, 5};
    tree.handle_event(ev1, sampler, rng);
    
    std::cout << "Tree after first insertion:\n" << tree.print_avl() << "\n";
    
    // Second insertion in the added part (position 51-54)
    Event ev2 = {INSERTION, 52, 3};
    print_event(ev2);
    
    bool success = tree.handle_event(ev2, sampler, rng);
    std::cout << "Event success: " << (success ? "YES" : "NO") << "\n";
    std::cout << "Tree after second insertion:\n" << tree.print_avl() << "\n";
    
    bool valid = tree.validate_rate_integrity();
    std::cout << "Rate integrity: " << (valid ? "PASS" : "FAIL") << "\n";
    
    return success && valid;
}

// Test 3: Multiple insertions
bool test_multiple_insertions() {
    print_test_header("Multiple Insertions");
    
    avl_array_with_rates<size_t, size_t, 1000> tree;
    std::vector<size_t> parentRates(100, 0);
    tree.init_tree(100, parentRates);
    
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    CategorySampler sampler(transMatrix, probs, MAX_INSERTION_LENGTH);
    std::mt19937_64 rng(42);
    
    std::vector<Event> events = {
        {INSERTION, 20, 3},
        {INSERTION, 50, 5},
        {INSERTION, 80, 2},
        {INSERTION, 25, 4}
    };
    
    bool all_valid = true;
    for (auto& ev : events) {
        print_event(ev);
        bool success = tree.handle_event(ev, sampler, rng);
        std::cout << "Event success: " << (success ? "YES" : "NO") << "\n";
        std::cout << "Tree state:\n" << tree.print_avl() << "\n";
        
        bool valid = tree.validate_rate_integrity();
        std::cout << "Rate integrity: " << (valid ? "PASS" : "FAIL") << "\n";
        
        all_valid = all_valid && success && valid;
    }
    
    return all_valid;
}

// Test 4: Deletion in original part (Case A)
bool test_deletion_case_a() {
    print_test_header("Deletion Case A - Delete from beginning of OP");
    
    avl_array_with_rates<size_t, size_t, 1000> tree;
    std::vector<size_t> parentRates(100, 0);
    tree.init_tree(100, parentRates);
    
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    CategorySampler sampler(transMatrix, probs, MAX_INSERTION_LENGTH);
    std::mt19937_64 rng(42);
    
    // Create a block with insertion
    Event ins = {INSERTION, 50, 10};
    tree.handle_event(ins, sampler, rng);
    std::cout << "After insertion:\n" << tree.print_avl() << "\n";
    
    // Delete from position 1 (first deletable position)
    Event del = {DELETION, 1, 20};  // Changed from 0 to 1
    print_event(del);
    bool success = tree.handle_event(del, sampler, rng);
    std::cout << "Event success: " << (success ? "YES" : "NO") << "\n";
    std::cout << "Tree after deletion:\n" << tree.print_avl() << "\n";
    
    bool valid = tree.validate_rate_integrity();
    std::cout << "Rate integrity: " << (valid ? "PASS" : "FAIL") << "\n";
    
    return success && valid;
}

// Test 5: Deletion spanning OP and AP (Case C)
bool test_deletion_case_c() {
    print_test_header("Deletion Case C - Delete spanning OP and AP");
    
    avl_array_with_rates<size_t, size_t, 1000> tree;
    std::vector<size_t> parentRates(100, 0);
    tree.init_tree(100, parentRates);
    
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    CategorySampler sampler(transMatrix, probs, MAX_INSERTION_LENGTH);
    std::mt19937_64 rng(42);
    
    // Create a block with insertion
    Event ins = {INSERTION, 50, 20};
    tree.handle_event(ins, sampler, rng);
    std::cout << "After insertion:\n" << tree.print_avl() << "\n";
    
    // Delete spanning OP and part of AP from position 1
    Event del = {DELETION, 1, 60}; // Deletes 50 from OP, 10 from AP
    print_event(del);
    bool success = tree.handle_event(del, sampler, rng);
    std::cout << "Event success: " << (success ? "YES" : "NO") << "\n";
    std::cout << "Tree after deletion:\n" << tree.print_avl() << "\n";
    
    bool valid = tree.validate_rate_integrity();
    std::cout << "Rate integrity: " << (valid ? "PASS" : "FAIL") << "\n";
    
    return success && valid;
}
// Test 6: Deletion in middle of AP (Case F)
bool test_deletion_case_f_in_ap() {
    print_test_header("Deletion Case F - Delete in middle of AP");
    
    avl_array_with_rates<size_t, size_t, 1000> tree;
    std::vector<size_t> parentRates(100, 0);
    tree.init_tree(100, parentRates);
    
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    CategorySampler sampler(transMatrix, probs, MAX_INSERTION_LENGTH);
    std::mt19937_64 rng(42);
    
    // Create a block with large insertion
    Event ins = {INSERTION, 50, 30};
    tree.handle_event(ins, sampler, rng);
    std::cout << "After insertion (block has length=50, insertion=30):\n" << tree.print_avl() << "\n";
    
    // Delete in the middle of the added part
    // Position 60 = 10 positions into the AP
    Event del = {DELETION, 60, 10};
    print_event(del);
    bool success = tree.handle_event(del, sampler, rng);
    std::cout << "Event success: " << (success ? "YES" : "NO") << "\n";
    std::cout << "Tree after deletion:\n" << tree.print_avl() << "\n";
    
    bool valid = tree.validate_rate_integrity();
    std::cout << "Rate integrity: " << (valid ? "PASS" : "FAIL") << "\n";
    
    return success && valid;
}

// Test 7: Complex sequence of insertions and deletions
bool test_complex_sequence() {
    print_test_header("Complex Sequence - Mixed Insertions and Deletions");
    
    avl_array_with_rates<size_t, size_t, 1000> tree;
    std::vector<size_t> parentRates(100, 0);
    tree.init_tree(100, parentRates);
    
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    CategorySampler sampler(transMatrix, probs, MAX_INSERTION_LENGTH);
    std::mt19937_64 rng(42);
    
    std::vector<Event> events = {
        {INSERTION, 20, 5},
        {INSERTION, 50, 10},
        {DELETION, 30, 15},
        {INSERTION, 40, 8},
        {DELETION, 10, 5},
        {INSERTION, 70, 3}
    };
    
    bool all_valid = true;
    for (size_t i = 0; i < events.size(); ++i) {
        auto& ev = events[i];
        print_event(ev);
        std::cout << "Step " << (i + 1) << "/" << events.size() << "\n";
        
        bool success = tree.handle_event(ev, sampler, rng);
        std::cout << "Event success: " << (success ? "YES" : "NO") << "\n";
        std::cout << "Tree state:\n" << tree.print_avl() << "\n";
        
        bool valid = tree.validate_rate_integrity();
        std::cout << "Rate integrity: " << (valid ? "PASS" : "FAIL") << "\n";
        
        if (!success || !valid) {
            std::cout << "ERROR: Test failed at step " << (i + 1) << "\n";
            all_valid = false;
            break;
        }
    }
    
    return all_valid;
}

int main() {
    std::cout << "AVL Tree with Rate Categories - Test Suite\n";
    std::cout << "==========================================\n";
    
    int passed = 0;
    int total = 7;
    
    if (test_single_insertion_in_original()) passed++;
    if (test_single_insertion_in_added()) passed++;
    if (test_multiple_insertions()) passed++;
    if (test_deletion_case_a()) passed++;
    if (test_deletion_case_c()) passed++;
    if (test_deletion_case_f_in_ap()) passed++;
    if (test_complex_sequence()) passed++;
    
    std::cout << "\n========================================\n";
    std::cout << "FINAL RESULTS: " << passed << "/" << total << " tests passed\n";
    std::cout << "========================================\n";
    
    return (passed == total) ? 0 : 1;
}