#ifndef CACHED_TRANSITION_PROBABILITIES_H
#define CACHED_TRANSITION_PROBABILITIES_H

#include <vector>
#include <unordered_map>
#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/DiscreteNDistribution.h"
#include "../libs/Phylolib/includes/stochasticProcess.h"

template<size_t AlphabetSize>
class CachedTransitionProbabilities {
public:
    CachedTransitionProbabilities(const tree& _tree, const stochasticProcess& _sp)
    {
        const size_t numNodes = _tree.getNodesNum();
        const size_t numCategories = _sp.categories();
        
        // Reserve mapping vector
        _nodeToUniqueIndex.resize(numNodes);
        
        // Map for deduplicating branch lengths
        std::unordered_map<long, size_t> branchLengthToIndex;
        
        // Iterative tree traversal
        std::vector<tree::nodeP> nodesToProcess;
        nodesToProcess.push_back(_tree.getRoot());
        
        size_t pos = 0;
        while (pos < nodesToProcess.size()) {
            tree::nodeP currentNode = nodesToProcess[pos];
            
            // Process non-root nodes
            if (pos > 0) {
                const MDOUBLE branchLength = currentNode->dis2father();
                const int nodeID = currentNode->id();
                const long key = static_cast<long>(branchLength * 1e6);
                
                // Check if this branch length already exists
                auto it = branchLengthToIndex.find(key);
                if (it == branchLengthToIndex.end()) {
                    // New unique branch length - create distributions
                    size_t newIndex = _distributions.size();
                    branchLengthToIndex[key] = newIndex;
                    
                    std::vector<DiscreteNDistribution<AlphabetSize>> nodeDistributions;
                    nodeDistributions.reserve(numCategories * AlphabetSize);
                    
                    for (size_t cat = 0; cat < numCategories; ++cat) {
                        const MDOUBLE rate = _sp.rates(cat);
                        
                        for (size_t i = 0; i < AlphabetSize; ++i) {
                            std::vector<double> probabilities;
                            probabilities.reserve(AlphabetSize);
                            MDOUBLE normalizingFactor = 0.0;
                            
                            for (size_t j = 0; j < AlphabetSize; ++j) {
                                MDOUBLE prob = _sp.Pij_t(i, j, branchLength * rate);
                                // std::cout << "(" << i << j << ")=" <<prob << "\n";

                                normalizingFactor += prob;
                                probabilities.push_back(prob);
                            }
                            
                            nodeDistributions.emplace_back((probabilities));
                        }
                    }
                    
                    _distributions.push_back(std::move(nodeDistributions));
                    _nodeToUniqueIndex[nodeID] = newIndex;
                } else {
                    // Reuse existing distributions
                    _nodeToUniqueIndex[nodeID] = it->second;
                }
            }
            
            // Add children to queue
            for (size_t k = 0; k < currentNode->getNumberOfSons(); ++k) {
                nodesToProcess.push_back(currentNode->getSon(k));
            }
            ++pos;
        }
    }
    
    DiscreteNDistribution<AlphabetSize>& getDistribution(int nodeID, int category, int character)  {
        size_t uniqueIndex = _nodeToUniqueIndex[nodeID];
        size_t distributionIndex = category * AlphabetSize + character;

        return _distributions[uniqueIndex][distributionIndex];
    }



    size_t getNumUniqueBranches() {return _distributions.size();}

private:
    // std::vector<std::vector<DiscreteNDistribution<20>>> _distributionsAmino;
    std::vector<std::vector<DiscreteNDistribution<AlphabetSize>>> _distributions;

    std::vector<size_t> _nodeToUniqueIndex;
};

#endif // CACHED_TRANSITION_PROBABILITIES_H