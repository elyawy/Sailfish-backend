#include <iostream>
#include <memory>
#include <cassert>
#include "../../../src/CachedTransitionProbabilities.h"
#include "../../../src/allModels.h"

#include "../../../libs/Phylolib/includes/tree.h"
#include "../../../libs/Phylolib/includes/stochasticProcess.h"
#include "../../../libs/Phylolib/includes/gammaDistribution.h"
#include "../../../libs/Phylolib/includes/nucJC.h"
#include "../../../libs/Phylolib/includes/trivialAccelerator.h"
#include "../../../libs/Phylolib/includes/aaJC.h"

int main() {
    // 1. Load tree from file
    tree testTree("../../trees/normalbranches_nLeaves10.treefile");
    std::cout << "Loaded tree with " << testTree.getNodesNum() << " nodes\n";
    
    // 2. Create stochastic process for amino acids with JC model
    auto repModel = std::make_unique<pupAll>(datMatrixHolder::wag);

    auto pij = std::make_unique<trivialAccelerator>(repModel.get());
    
    const MDOUBLE alpha = 1.0;
    const int gammaCategories = 4;
    gammaDistribution dist(alpha, gammaCategories);
    
    auto sp = std::make_shared<stochasticProcess>(&dist, pij.get());
    
    std::cout << "Created stochastic process with " 
              << sp->categories() << " categories and alphabet size " 
              << sp->alphabetSize() << "\n";
    
    // 3. Create cached transition probabilities
    CachedTransitionProbabilities<20> cachedPij(testTree, *sp);
    std::cout << "Built CachedTransitionProbabilities\n";
    
    // 4. Test retrieval - get distribution for node 1, category 0, character 0
    if (testTree.getNodesNum() > 1) {
        const DiscreteNDistribution<20>& dist = cachedPij.getDistribution(3, 2, 5);
        dist.printTable();
        std::cout << "Successfully retrieved distribution for node 3 category 2 and character 5\n";
    }
    
    std::cout << "Test completed successfully!\n";
    return 0;
}