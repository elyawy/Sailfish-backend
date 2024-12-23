#include <iostream>
#include "../../src/BlockTree.h"
#include "../../src/Sequence.h"
#include "../../src/MSA.h"

int main() {

    size_t root_size = 5;
    size_t random_position = root_size+1;
    SuperSequence superSeq(root_size, 5);
    superSeq.printSequence();
    
    Sequence origin_seq(superSeq, 1, 0);


    origin_seq.initSequence();
    origin_seq.printSequence();

    BlockTree tree(root_size);
    tree.handleEvent(event::INSERTION, 2,2);
    tree.handleEvent(event::INSERTION, 7,3);

    std::cout << tree.printTree() << "\n";

    Sequence seq(superSeq, 1, 1);
    seq.generateSequence(tree.getBlockList(), origin_seq);
    superSeq.printSequence();
    seq.printSequence();


    BlockTree tree2(10);
    tree2.handleEvent(event::INSERTION, 1,4);
    tree2.handleEvent(event::INSERTION, 5,2);

    std::cout << tree2.printTree() << "\n";

    Sequence seq2(superSeq, 1 , 2);
    seq2.generateSequence(tree2.getBlockList(), seq);
    superSeq.printSequence();
    seq2.printSequence();


    BlockTree tree3(10);
    tree3.handleEvent(event::INSERTION, 4,1);
    tree3.handleEvent(event::DELETION, 5,3);
    std::cout << tree3.printTree() << "\n";

    Sequence seq3(superSeq, 1, 3);
    seq3.generateSequence(tree3.getBlockList(), seq);
    superSeq.printSequence();
    seq3.printSequence();

    std::vector<Sequence> finalSeqs;
    
    finalSeqs.push_back(origin_seq);
    finalSeqs.push_back(seq);
    finalSeqs.push_back(seq2);
    finalSeqs.push_back(seq3);

    BlockTree righttree(root_size);
    righttree.handleEvent(event::DELETION, 2,1);
    righttree.handleEvent(event::INSERTION, 5,1);
    righttree.handleEvent(event::INSERTION, 4,3);
    std::cout << righttree.printTree() << "\n";


    Sequence rightseq(superSeq, 1, 4);
    rightseq.generateSequence(righttree.getBlockList(), origin_seq);
    superSeq.printSequence();
    rightseq.printSequence();

    finalSeqs.push_back(rightseq);


    MSA msa(finalSeqs.size(), superSeq.getMsaSequenceLength(), {1,1,1,1,1});
    msa.fillMSA(finalSeqs, superSeq);
    std::cout << msa.generateMsaStringWithoutSubs();

    // for(auto item: superSequence) {
    //     std::cout << item << " ";
    // }
    // std::cout << "\n";


    return 0;
}