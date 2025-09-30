#include <iostream>
#include "../../src/LinkedList.h"


int main() {

    FixedList fl(20);

    fl.initialize(10);

    size_t index_1 = fl.insertAfter(4, true);
    size_t index_2 = fl.insertAfter(10, true);

    fl.referencePosition(6);

    std::cout << index_1 << " " << index_2 << "\n";

    fl.setAbsolutePositions();

    fl.printSequence();

    fl.printIndicesVector();

    fl.printTraversalVec();

    fl.checkSequenceValidity();



    std::cout << fl.getAbsolutePosition(index_1) << " ";
    std::cout << fl.getAbsolutePosition(index_2) << " ";

    std::cout << "\n";



    return 0;
}