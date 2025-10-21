#include "../../src/Sequence.h"
#include "../../src/FixedList.h"
#include "../../src/BlockTree.h"
#include <iostream>
#include <cassert>

void printSequence(const Sequence& seq, const std::string& name) {
    std::cout << name << " positions: ";
    seq.printSequence();
}

// Scenario 1: Linear chain - Root -> Child1 -> Child2
void testLinearChain() {
    std::cout << "=== Scenario 1: Linear Chain (Root -> Child1 -> Child2) ===\n\n";
    BlockTree bt;
    size_t root_size = 10;
    
    std::cout << "create super seq\n";

    // Create super sequence with space for 50 elements
    FixedList superSeq(50);
    superSeq.initialize(root_size);
    
    std::cout << "generate root\n";

    // Generate root sequence
    Sequence rootSeq(superSeq, true, 0);
    std::cout << "init root\n";

    rootSeq.initSequence();
    
    std::cout << "Root sequence initialized with " << root_size << " positions\n";
    printSequence(rootSeq, "Root");
    if (!rootSeq.checkSequenceValidity()) std::cout << "Root invalid!\n";
    
    superSeq.setAbsolutePositions();
    std::cout << "Root MSA length: " << superSeq.getMsaSequenceLength() << "\n\n";
    
    // Generate child1 from root: Delete 3 positions starting at position 5
    bt.initTree(root_size);
    // BlockTree btree1(root_size);
    bt.handleEvent(event::DELETION, 5, 3);
    std::cout << "Child1 BlockTree:\n" << bt.printTree() << "\n";
    
    Sequence child1Seq(superSeq, true, 1);
    child1Seq.generateSequence(bt.getBlockList(), rootSeq);
    
    std::cout << "Child1: deleted 3 positions starting at position 5\n";
    printSequence(child1Seq, "Child1");
    if (!child1Seq.checkSequenceValidity()) std::cout << "Child1 invalid!\n";
    
    superSeq.setAbsolutePositions();
    std::cout << "After Child1 MSA length: " << superSeq.getMsaSequenceLength() << "\n\n";
    
    bt.clear();
    // Generate child2 from child1: Insert 2 positions at position 3
    bt.initTree(child1Seq.size());
    bt.handleEvent(event::INSERTION, 3, 2);
    std::cout << "Child2 BlockTree:\n" << bt.printTree() << "\n";
    
    Sequence child2Seq(superSeq, true, 2);
    child2Seq.generateSequence(bt.getBlockList(), child1Seq);
    
    std::cout << "Child2: inserted 2 positions at position 3\n";
    printSequence(child2Seq, "Child2");
    if (!child2Seq.checkSequenceValidity()) std::cout << "Child2 invalid!\n";
    
    superSeq.setAbsolutePositions();
    std::cout << "After Child2 MSA length: " << superSeq.getMsaSequenceLength() << "\n\n";
    
    // Verify sequence lengths
    std::cout << "Root size: " << rootSeq.size() << "\n";
    std::cout << "Child1 size: " << child1Seq.size() << "\n";
    std::cout << "Child2 size: " << child2Seq.size() << "\n";
    
    assert(rootSeq.size() == root_size);  // +1 for anchor
    assert(child1Seq.size() == rootSeq.size() - 3);  // -3 deletions
    assert(child2Seq.size() == child1Seq.size() + 2);  // +2 insertions
    
    std::cout << "Scenario 1 passed!\n\n";
}

// Scenario 2: Branching - Root with two different children
void testBranching() {
    std::cout << "=== Scenario 2: Branching (Root -> Child1 and Child2) ===\n\n";
    BlockTree bt;
    size_t root_size = 10;
    
    std::cout << "create super seq\n";
    
    // Create super sequence with space for 50 elements
    FixedList superSeq(50);
    superSeq.initialize(root_size);
    
    std::cout << "generate root\n";
    
    // Generate root sequence
    Sequence rootSeq(superSeq, true, 0);
    std::cout << "init root\n";
    
    rootSeq.initSequence();
    
    std::cout << "Root sequence initialized with " << root_size << " positions\n";
    printSequence(rootSeq, "Root");
    if (!rootSeq.checkSequenceValidity()) std::cout << "Root invalid!\n";
    
    superSeq.setAbsolutePositions();
    std::cout << "Root MSA length: " << superSeq.getMsaSequenceLength() << "\n\n";
    
    // Generate child1 from root: Insert 3 positions at position 4
    bt.initTree(root_size);
    bt.handleEvent(event::INSERTION, 4, 3);
    std::cout << "Child1 BlockTree:\n" << bt.printTree() << "\n";
    
    Sequence child1Seq(superSeq, true, 1);
    child1Seq.generateSequence(bt.getBlockList(), rootSeq);
    
    std::cout << "Child1: inserted 3 positions at position 4\n";
    printSequence(child1Seq, "Child1");
    if (!child1Seq.checkSequenceValidity()) std::cout << "Child1 invalid!\n";
    
    superSeq.setAbsolutePositions();
    std::cout << "After Child1 MSA length: " << superSeq.getMsaSequenceLength() << "\n\n";
    
    bt.clear();
    
    // Generate child2 from root: Delete 4 positions starting at position 2
    bt.initTree(root_size);
    bt.handleEvent(event::DELETION, 2, 4);
    std::cout << "Child2 BlockTree:\n" << bt.printTree() << "\n";
    
    Sequence child2Seq(superSeq, true, 2);
    child2Seq.generateSequence(bt.getBlockList(), rootSeq);
    
    std::cout << "Child2: deleted 4 positions starting at position 2\n";
    printSequence(child2Seq, "Child2");
    if (!child2Seq.checkSequenceValidity()) std::cout << "Child2 invalid!\n";
    
    superSeq.setAbsolutePositions();
    std::cout << "After Child2 MSA length: " << superSeq.getMsaSequenceLength() << "\n\n";
    
    // Verify sequence lengths
    std::cout << "Root size: " << rootSeq.size() << "\n";
    std::cout << "Child1 size: " << child1Seq.size() << "\n";
    std::cout << "Child2 size: " << child2Seq.size() << "\n";
    
    assert(rootSeq.size() == root_size);
    assert(child1Seq.size() == rootSeq.size() + 3);  // +3 insertions
    assert(child2Seq.size() == rootSeq.size() - 4);  // -4 deletions
    
    std::cout << "Scenario 2 passed!\n\n";
}

int main() {
    testLinearChain();
    testBranching();
    
    std::cout << "All tests completed!\n";
    return 0;
}