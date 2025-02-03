#include <iostream>
#include "../../src/BlockTree.h"
#include "../../src/Sequence.h"
#include "../../src/MSA.h"

int main() {

    size_t root_size = 20;
    size_t random_position = root_size+1;
    SuperSequence superSeq(root_size, 1);
    // superSeq.printSequence();
    
    Sequence origin_seq(superSeq, true, 0);
    origin_seq.initSequence();
    origin_seq.printSequence();

    BlockTree tree(root_size);
    // tree.handleEvent(event::DELETION, 0, 100);
    tree.handleEvent(event::INSERTION, 0, 10);

    std::cout << tree.printTree() << "\n";

    Sequence seq(superSeq, true, 1);
    seq.generateSequence(tree.getBlockList(), origin_seq);
    // superSeq.printSequence();
    seq.printSequence();
    if (!seq.checkSequenceValidity()) std::cout << "invalid\n";

    // for(auto item: superSequence) {
    //     std::cout << item << " ";
    // }
    // std::cout << "\n";


    return 0;
}