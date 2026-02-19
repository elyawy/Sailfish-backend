#include <iostream>
#include <cassert>

#include "../../../src/AvlTree.h"

// Helper to build the evolved sequence string from the tree for verification.
std::string build_sequence(avl_array<size_t, size_t, 100> &tree) {
    std::string seq = "";
    size_t last_position = 0;
    for (auto it = tree.begin(); it != tree.end(); ++it) {
        seq += std::string(it.key() - last_position, 'X');
        if (it.key() != 0) {
            seq += std::string(it.val().length, 'R');
        } else {
            seq += std::string(1, 'A');
            seq += std::string(it.val().length - 1, 'R');
        }
        seq += std::string(it.val().insertion, 'N');
        last_position = it.key() + it.val().length; // Adjust last_position to account for insertions
    }
    return seq;
}

int main() {
    // Starting sequence of length 10: A R R R R R R R R R
    // Positions:                       0 1 2 3 4 5 6 7 8 9
    avl_array<size_t, size_t, 100> tree;
    tree.init_tree(10);

    // --- Event 1: Insertion at position 3 (size 2) ---
    // Inserts 2 new bases after position 3 (inside the original block).
    // Expected: A R R R N N R R R R R R
    size_t pos1 = 3;
    tree.handle_event(INSERTION, pos1, 2);
    std::cout << "[After insertion at 3, size 2]\n";
    std::cout << "Tree:\n" << tree.print_avl();
    std::string seq1 = build_sequence(tree);
    std::cout << "Sequence: " << seq1 << "\n\n";
    assert(tree.checkLength());
    assert(seq1 == "ARRRNNRRRRRR");

    // --- Event 2: Deletion at position 0 (size 1) ---
    // Deletes 1 base from the very start of the sequence.
    // Position 0 is the anchor 'A', which is preserved as a stub.
    // Expected: A X R R N N R R R R R R
    size_t pos2 = 1;
    tree.handle_event(DELETION, pos2, 1);
    std::cout << "[After deletion at 1, size 1]\n";
    std::cout << "Tree:\n" << tree.print_avl();
    std::string seq2 = build_sequence(tree);
    std::cout << "Sequence: " << seq2 << "\n\n";
    assert(tree.checkLength());
    assert(seq2 == "AXRRNNRRRRRR");

    // --- Event 3: Deletion at position 4 (size 3) ---
    // Deletes 3 bases starting at position 4, spanning into the insertion region.
    // This deletion should span the boundary between the two blocks.
    // Expected: A X R R N X X R R R R
    size_t pos3 = 4;
    tree.handle_event(DELETION, pos3, 3);
    std::cout << "[After deletion at 4, size 3]\n";
    std::cout << "Tree:\n" << tree.print_avl();
    std::string seq3 = build_sequence(tree);
    std::cout << "Sequence: " << seq3 << "\n\n";
    assert(tree.checkLength());
    assert(seq3 == "AXRRNXXRRRR");

    // --- Event 4: Insertion at position 6 (size 2) ---
    // Inserts 2 bases near the right edge of the original sequence.
    // Expected: A X R R N X X R R N N R
    size_t pos4 = 6;
    tree.handle_event(INSERTION, pos4, 2);
    std::cout << "[After insertion at 6, size 2]\n";
    std::cout << "Tree:\n" << tree.print_avl();
    std::string seq4 = build_sequence(tree);
    std::cout << "Sequence: " << seq4 << "\n\n";
    assert(tree.checkLength());
    assert(seq4 == "AXRRNXXRRRNNR");

    // --- Event 5: Deletion at position 7 (size 4) ---
    // Deletes 4 bases starting at position 7, spanning into the insertion block
    // and the following original block.
    // Expected: A X R R X X X R X X X X
    size_t pos5 = 7;
    tree.handle_event(DELETION, pos5, 4);
    std::cout << "[After deletion at 7, size 4]\n";
    std::cout << "Tree:\n" << tree.print_avl();
    std::string seq5 = build_sequence(tree);
    std::cout << "Sequence: " << seq5 << "\n\n";
    assert(tree.checkLength());
    assert(seq5 == "AXRRNXXRRR");

    std::cout << "All assertions passed!\n";
    return 0;
}