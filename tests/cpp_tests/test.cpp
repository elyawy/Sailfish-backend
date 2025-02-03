#include <iostream>
#include "../../src/BlockTree.h"
#include "../../src/Sequence.h"
#include "../../src/MSA.h"

int main() {

    size_t root_size = 10;
    size_t random_position = root_size+1;
    SuperSequence superSeq(root_size, 3);
    // superSeq.printSequence();

    
    Sequence origin_seq(superSeq, true, 0);
    origin_seq.initSequence();
    origin_seq.printSequence();

    BlockTree btree(root_size);
    btree.handleEvent(event::DELETION, 1, 10);
    std::cout << btree.printTree() << "\n";
    Sequence seq(superSeq, true, 1);
    seq.generateSequence(btree.getBlockList(), origin_seq);
    seq.printSequence();
    if (!seq.checkSequenceValidity()) std::cout << "invalid\n";

    BlockTree btree2(0);
    btree2.handleEvent(event::INSERTION, 0, 1);
    std::cout << btree2.printTree() << "\n";
    Sequence seq2(superSeq, true, 2);
    seq2.generateSequence(btree2.getBlockList(), seq);
    seq2.printSequence();
    if (!seq2.checkSequenceValidity()) std::cout << "invalid\n";


    std::vector<Sequence> seqs;

    seqs.push_back(origin_seq);
    seqs.push_back(seq);
    seqs.push_back(seq2);


    MSA msa = MSA::msaFromSequences(seqs, superSeq);
    std::cout << msa.generateMsaStringWithoutSubs() << "\n";

    
    // for(auto item: superSequence) {
    //     std::cout << item << " ";
    // }
    // std::cout << "\n";


    return 0;
}