#include <iostream>
#include <cassert>

#include "../../../src/AvlTree.h"

int main() {
    // simple test to check the way events are mapped to the tree structure.
    avl_array<size_t, size_t, 100> tree;
    tree.init_tree(10);

    size_t event_position_1 = 0;
    tree.handle_event(DELETION, event_position_1, 1);

    size_t event_position_2 = 5; 
    tree.handle_event(DELETION, event_position_2, 1);

    size_t event_position_3 = 5; 
    tree.handle_event(DELETION, event_position_3, 1);

    size_t event_position_4 = 6; 
    tree.handle_event(DELETION, event_position_4, 1);

    size_t event_position_5 = 2; 
    tree.handle_event(DELETION, event_position_5, 1);

    std::string evolved_sequence = "";
    size_t last_position = 0;
    for (auto it = tree.begin(); it != tree.end(); ++it) {
        std::cout << "Position: " << it.key() << ", Length: " << it.val().length << ", Insertion: " << it.val().insertion << "\n";
        evolved_sequence += std::string(it.key() - last_position, 'X');
        if (it.key() != 0) {
            evolved_sequence += std::string(it.val().length, 'R');
        } else {
            evolved_sequence += std::string(1, 'A');
            evolved_sequence += std::string(it.val().length-1, 'R');
        }
        evolved_sequence += std::string(it.val().insertion, 'N');
        last_position = it.key() + it.val().length;
    }

    assert(evolved_sequence == "ARXRRXXRXR");

    std::cout << "Evolved sequence: " << evolved_sequence << "\n";

    return 0;
}