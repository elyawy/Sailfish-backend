#ifndef ___SIMULATION_CONTEXT_H
#define ___SIMULATION_CONTEXT_H

#include <vector>
#include <random>
#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/sequence.h"
#include "SimulationProtocol.h"
#include "CategorySampler.h"

typedef std::string SparseSequence;
typedef std::vector<SparseSequence> SparseSequenceContainer;



constexpr uint64_t PHI = 0x9e3779b97f4a7c15;

template<typename RngType = std::mt19937_64>
class SimulationContext {
public:
    SimulationContext(tree* t, size_t seed, SimulationProtocol* protocol = nullptr) 
        : _tree(t), _seed(seed), _rng(seed*PHI), 
          _nodesToSave(t->getNodesNum(), false), 
          _idToSaveIndices(t->getNodesNum(), SIZE_MAX),
          _protocol(protocol) {
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


    
    void setSaveRoot() { 
        if (_nodesToSave[_tree->getRoot()->id()]) return; // already set
        _nodesToSave[_tree->getRoot()->id()] = true;
        // since root is the first node it should be first in the indices and names lists
        _nodesToSaveIndices.insert(_nodesToSaveIndices.begin(), _tree->getRoot()->id());
        _idToSaveIndices[_tree->getRoot()->id()] = 0;
        _nodeToSaveNames.insert(_nodeToSaveNames.begin(), _tree->getRoot()->name());
        _numberOfNodesToSave++;
    }

    void setSaveLeaves() {
        // fill bools, indices and names
        _numberOfNodesToSave = 0;
        _nodesToSaveIndices.clear();
        _nodeToSaveNames.clear();
        setSaveLeavesRecursive(_tree->getRoot());
    }

    void setSaveAll() {
        // fill bools, indices and names
        _numberOfNodesToSave = 0;
        _nodesToSaveIndices.clear();
        _nodeToSaveNames.clear();
        setAllNodesRecursive(_tree->getRoot());
    }

    const std::vector<bool>& getNodesToSave() const { return _nodesToSave; }
    const std::vector<size_t>& getIdToSaveIndices() const { return _idToSaveIndices; }

    size_t getNumberOfNodesToSave() const { return _numberOfNodesToSave; }

    const std::vector<size_t>& getNodesToSaveIndices() const { return _nodesToSaveIndices; }
    const std::vector<std::string>& getNodeToSaveNames() const { return _nodeToSaveNames; }

    const tree::nodeP getRootNode() const { return _tree->getRoot(); }

    void setProtocol(SimulationProtocol* protocol) {
        _protocol = protocol;
    }

        // Getter - returns nullptr if not set
    SimulationProtocol* getProtocol() const { 
        return _protocol; 
    }

    void setCategorySampler(CategorySampler* categorySampler) {
        _categorySampler = categorySampler;
    }

    CategorySampler* getCategorySampler() {
        return _categorySampler;
    }

private:
    void setSaveLeavesRecursive(tree::nodeP node) {
        if (node->isLeaf()) {
            setNodeToSaveInfo(node);
            return;
        }
        // Not a leaf - do not save
        _nodesToSave[node->id()] = false;
        _idToSaveIndices[node->id()] = SIZE_MAX;
        for (auto& child : node->getSons()) {
            setSaveLeavesRecursive(child);
        }
    }

    void setAllNodesRecursive(tree::nodeP node) {
        setNodeToSaveInfo(node);
        for (auto& child : node->getSons()) {
            setAllNodesRecursive(child);
        }
    }

    void setNodeToSaveInfo(tree::nodeP node) {
        _nodesToSave[node->id()] = true;
        _idToSaveIndices[node->id()] = _numberOfNodesToSave;
        _nodesToSaveIndices.push_back(node->id());
        _nodeToSaveNames.push_back(node->name());
        _numberOfNodesToSave++;
    }



    tree* _tree;
    size_t _seed;
    RngType _rng;
    std::vector<bool> _nodesToSave;
    std::vector<size_t> _nodesToSaveIndices;
    std::vector<size_t> _idToSaveIndices; // maps from node id to index in saved nodes. if a node id is not saved, its index is SIZE_MAX
    std::vector<std::string> _nodeToSaveNames;
    size_t _numberOfNodesToSave = 0;
    SimulationProtocol* _protocol;  // Can be nullptr
    CategorySampler* _categorySampler; // Can be nullptr
};

#endif