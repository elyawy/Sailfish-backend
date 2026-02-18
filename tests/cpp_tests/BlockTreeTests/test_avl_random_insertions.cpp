#include <iostream>
#include <vector>
#include <random>
#include "../../../src/BlockCommon.h"
#include "../../../src/AvlTreeWithRates.h"
#include "../../../src/Event.h"
#include "../../../src/CategorySampler.h"


const size_t MAX_INSERTION_LENGTH = 100;

CategorySampler create_sampler() {
    std::vector<MDOUBLE> probs = {0.25, 0.25, 0.25, 0.25};
    std::vector<std::vector<MDOUBLE>> transMatrix(4, probs);
    return CategorySampler(transMatrix, probs, MAX_INSERTION_LENGTH);
} 

// ten thousand random insertion events are sampled, between 0 to currentSequenceLength, with insertion length between 1 and 10 (uniform). After each insertion, the integrity of the tree is checked.
int main(){ 
    avl_array_with_rates<size_t, size_t, 100> tree;

    size_t initialSequenceLength = 50;

    std::vector<size_t> parentRates(initialSequenceLength, 0);
    tree.init_tree(initialSequenceLength + 1, parentRates);
    
    CategorySampler sampler = create_sampler();
    std::mt19937_64 rng(10);


    size_t currestSequenceLength = initialSequenceLength;
    const int NUM_INSERTIONS = 10000;
    for (int i = 0; i < NUM_INSERTIONS; i++) {
        std::uniform_int_distribution<size_t> posDist(0, currestSequenceLength);
        std::uniform_int_distribution<size_t> lenDist(1, MAX_INSERTION_LENGTH);
        
        size_t pos = posDist(rng);
        size_t len = lenDist(rng);
        
        Event ev = {INSERTION, pos, len};
        std::cout << "Insertion " << i+1 << ": position=" << pos << ", length=" << len << "\n";
        if (!tree.handle_event(ev, sampler, rng)) {
            std::cout << "FAIL: Insertion " << i+1 << " returned false\n";
            return false;
        }
        currestSequenceLength += len;
        
        if (!tree.validate_rate_integrity()) {
            std::cout << "FAIL: Rate integrity check failed after insertion " << i+1 << "\n";
            return false;
        }
    }
    
    std::cout << "PASS: All random insertions succeeded with valid tree integrity\n";
    return true;
}