#ifndef ___SIMULATION_CONTEXT_H
#define ___SIMULATION_CONTEXT_H

#include <vector>
#include <random>
#include "../libs/Phylolib/includes/tree.h"

constexpr uint64_t PHI = 0x9e3779b97f4a7c15;

template<typename RngType = std::mt19937_64>
class SimulationContext {
public:
    SimulationContext(tree* t, size_t seed) 
        : _tree(t), _seed(seed), _rng(seed*PHI), 
          _nodesToSave(t->getNodesNum(), false) {
            // By default, save leaves only
            setSaveLeaves();
          }

    tree* getTree() const { return _tree; }

    tree::nodeP getRoot() const { return _tree->getRoot(); }

    RngType& getRng() { return _rng; }
    
    size_t getSeed() const { return _seed; }
    
    void reseed(size_t seed) { 
        _seed = seed;
        _rng.seed(seed*PHI); 
    }

    const std::vector<bool>& getNodesToSave() const { return _nodesToSave; }

    void setSaveNode(size_t nodeId) { _nodesToSave[nodeId] = true; }
    
    void setSaveRoot() { _nodesToSave[_tree->getRoot()->id()] = true; }

    void setSaveLeaves() {
        setSaveLeavesRecursive(_tree->getRoot());
    }

    void setSaveAll() {
        std::fill(_nodesToSave.begin(), _nodesToSave.end(), true);
    }

private:
    void setSaveLeavesRecursive(tree::nodeP node) {
        if (node->isLeaf()) {
            _nodesToSave[node->id()] = true;
            return;
        }
        // Not a leaf - do not save
        _nodesToSave[node->id()] = false;
        for (auto& child : node->getSons()) {
            setSaveLeavesRecursive(child);
        }
    }

    tree* _tree;
    size_t _seed;
    RngType _rng;
    std::vector<bool> _nodesToSave;
};

#endif