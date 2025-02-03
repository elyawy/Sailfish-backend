#include <iostream>
#include "../../src/BlockTree.h"
#include "../../src/Sequence.h"
#include "../../src/MSA.h"

int main() {

    size_t root_size = 100;
    size_t random_position = root_size+1;
    SuperSequence superSeq(root_size, 2);
    // superSeq.printSequence();

    
    Sequence origin_seq(superSeq, true, 0);
    origin_seq.initSequence();
    origin_seq.printSequence();

    BlockTree btree(root_size);
    btree.handleEvent(event::DELETION, 1, 100);
    // tree.handleEvent(event::INSERTION, 11, 1);
    // tree.handleEvent(event::INSERTION, 12, 1);

    std::cout << btree.printTree() << "\n";

    Sequence seq(superSeq, true, 1);
    seq.generateSequence(btree.getBlockList(), origin_seq);
    // superSeq.printSequence();
    seq.printSequence();
    if (!seq.checkSequenceValidity()) std::cout << "invalid\n";

    std::vector<Sequence> seqs;

    seqs.push_back(origin_seq);
    seqs.push_back(seq);


    MSA msa = MSA::msaFromSequences(seqs, superSeq);
    std::cout << msa.generateMsaStringWithoutSubs() << "\n";

    
    // for(auto item: superSequence) {
    //     std::cout << item << " ";
    // }
    // std::cout << "\n";


    return 0;
}