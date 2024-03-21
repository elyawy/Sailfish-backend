
#include <memory>
#include <unordered_map>
#include <map>
#include <iostream>
#include <fstream>

#include "../libs/Phylolib/includes/sequence.h"



class substitutionManager
{
private:
    using changeMap = std::unordered_map<size_t, int>;
    std::vector<std::unique_ptr<changeMap>> _substitutionVec;
    // size_t _changeCounter;
public:
    substitutionManager(int numberOfTreeNodes) {
        _substitutionVec.resize(numberOfTreeNodes);

        
        // _changeCounter = 0;
    }

    int getCharacter(const int nodeId, const size_t position, const sequence &rootSeq) {
        if (_substitutionVec[nodeId] == nullptr) return rootSeq[position];
        changeMap* currentChanges = (_substitutionVec[nodeId].get());
        size_t isChanged = (*currentChanges).count(position);
        int returnChar =  isChanged ? (*currentChanges)[position] : rootSeq[position];
        return returnChar;
    }

    

    void handleEvent(const int nodeId, const size_t position, const int change) {
        // std::cout << "Change in node: " << nodeId << "\n";
        // std::cout << position << "->" << change << "\n";
        if (_substitutionVec[nodeId] == nullptr) {
            _substitutionVec[nodeId] = std::make_unique<changeMap>();
        }
        (*_substitutionVec[nodeId])[position] = change;
        
        // (_substitutionVec[nodeId])[position] = change;
        // ++_changeCounter;
        // if (_substitutionMap.count(nodeId) == 0) {
        //     std::shared_ptr<changeMap> seqChanges  = std::make_shared<changeMap>();
        //     (*seqChanges)[position] = change; // log change in changeMap;
        //     _substitutionMap[nodeId] = seqChanges;
        // } else {
        //     std::shared_ptr<changeMap> seqChanges = _substitutionMap[nodeId];
        //     (*seqChanges)[position] = change;
        //     _substitutionMap[nodeId] = seqChanges;
        // }
    }

    void dumpSubstitutionLog(const int nodeId) {
        if (nodeId == 0) return;
        if (_substitutionVec[nodeId] == nullptr) return;
        std::stringstream ss;
        std::stringstream changeString;

        ss << "/home/elyawy/temp/" << nodeId << ".sfasta";
        std::ofstream o(ss.str());

        for (auto &changes: *_substitutionVec[nodeId]) {
            changeString << changes.first << "-" << changes.second << "\n";
        }
        o << changeString.str();

    }

    // void readSubsLog(const int nodeId) {

    // }

    

    int getChange(const int nodeId, const size_t position) {
        if ((_substitutionVec[nodeId] == nullptr)) return -1;

        return (*_substitutionVec[nodeId])[position];
    }

    bool isEmpty(const int nodeId) {
        return _substitutionVec[nodeId] == nullptr;
    }


    std::unique_ptr<changeMap> getChangeMap(const int nodeId) {
        if ((_substitutionVec[nodeId] == nullptr)) {
            std::cout << nodeId << " <- this node is null\n";
        }

        return std::move(_substitutionVec[nodeId]);
    }

    void printSubManager() {
        std::cout << "printing subs...\n";
        for(auto &node: _substitutionVec) {
            for(auto &changes: *node){
                std::cout << changes.first << "->" << changes.second << ", ";
            }
            std::cout << "\n";
        }
    }

    // void printNumChangesResetCounter() {
        // std::cout << _changeCounter << "\n";
        // _changeCounter = 0;
    // };


    ~substitutionManager() {
        // for (size_t i = 0; i < _substitutionVec.size(); i++) {
        //     if (_substitutionVec[i] != nullptr)
        //     {
        //         delete _substitutionVec[i];
        //     }
        // }

    };
};
