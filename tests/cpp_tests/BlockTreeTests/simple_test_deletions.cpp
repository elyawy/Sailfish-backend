#include <iostream>


#include "../../../src/AvlTree.h"

int main() {
    // simple test to check the way events are mapped to the tree structure.
    avl_array<size_t, size_t, 100> tree;
    tree.init_tree(10);

    size_t event_position_1 = 1;

    std::cout << tree.print_avl();
    tree.handle_event(DELETION, event_position_1, 1);
    std::cout << "After deletion at position " << event_position_1 << " relative to parent" << ":\n";
    std::cout <<  tree.print_avl() << "\n";

    size_t event_position_2 = 1; 

    tree.handle_event(DELETION, event_position_2, 7);
    std::cout << "After deletion at position " << event_position_2 << ":\n";
    std::cout << tree.print_avl() << "\n";

    return 0;
}