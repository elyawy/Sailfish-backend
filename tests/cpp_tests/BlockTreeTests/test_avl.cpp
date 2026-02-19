
#include <iostream>
#include <vector>
#include <random>
#include "../../../src/AvlTree.h"
#include "../../../src/Event.h"





int main() {
    avl_array<size_t, size_t, 1000> tree;
    
    tree.init_tree(10);
    
    tree.handle_event(event::DELETION, 5, 1);
    std::cout << tree.print_avl();

    tree.handle_event(event::INSERTION, 6, 1);
    std::cout << tree.print_avl();

    tree.handle_event(event::INSERTION, 7, 1);
    std::cout << tree.print_avl();



    return 0;
}