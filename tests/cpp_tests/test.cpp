#include <iostream>
#include "../src/BlockTree.h"
#include "../src/Sequence.h"
#include "../src/MSA.h"

int main() {

    size_t root_size = 5;
    size_t random_position = root_size+1;
    SuperSequence superSeq(root_size);
    superSeq.printSequence();
    
    Sequence origin_seq(superSeq);


    origin_seq.initSequence();
    origin_seq.printSequence();

    BlockTree tree(root_size);
    tree.handleEvent(event::INSERTION, 2,5);
    tree.handleEvent(event::INSERTION, 9,3);


    Sequence seq(superSeq);
    seq.generateSequence(tree, origin_seq);
    superSeq.printSequence();
    seq.printSequence();

    BlockTree tree2(13);
    tree2.handleEvent(event::INSERTION, 1,4);
    tree2.handleEvent(event::INSERTION, 5,2);

    Sequence seq2(superSeq);
    seq2.generateSequence(tree2, seq);
    superSeq.printSequence();
    seq2.printSequence();


    BlockTree tree3(19);
    tree3.handleEvent(event::DELETION, 5,10);
    tree3.handleEvent(event::INSERTION, 1,6);

    Sequence seq3(superSeq);
    seq3.generateSequence(tree3, seq2);
    superSeq.printSequence();
    seq3.printSequence();

    std::vector<Sequence> finalSeqs;
    
    finalSeqs.push_back(origin_seq);
    finalSeqs.push_back(seq);
    finalSeqs.push_back(seq2);
    finalSeqs.push_back(seq3);

    MSA msa(finalSeqs);

    msa.printMSA();

    // for(auto item: superSequence) {
    //     std::cout << item << " ";
    // }
    // std::cout << "\n";


    return 0;
}