
#include <memory>
#include <unordered_map>
#include <map>
#include <iostream>
#include <fstream>

#include "../libs/Phylolib/includes/sequence.h"


const unsigned char INVALID_CHAR = 255;


class substitutionManager
{
private:
    using changeMap = std::vector<ALPHACHAR>;
    std::vector<std::unique_ptr<changeMap>> _substitutionVec;
    MDOUBLE _sumOfReactantsXRates;
    // size_t _changeCounter;
public:
    substitutionManager(int numberOfTreeNodes) {
        _substitutionVec.resize(numberOfTreeNodes);
        _sumOfReactantsXRates = 0.0;
        // _changeCounter = 0;
    }

    void updateReactantsSum(MDOUBLE Qii, MDOUBLE siteRate) {
        _sumOfReactantsXRates += (-Qii*siteRate);
        if (_sumOfReactantsXRates < 0) _sumOfReactantsXRates = 0.0;
    }

    MDOUBLE getReactantsSum() {
        return _sumOfReactantsXRates;
    }


    ALPHACHAR getCharacter(const int nodeId, const size_t position, const sequence &rootSeq) {
        if (_substitutionVec[nodeId] == nullptr) return rootSeq[position];
        

        changeMap* currentChanges = (_substitutionVec[nodeId].get());
        // size_t isChanged = (*currentChanges).count(position);
        // size_t isChanged = ((*currentChanges)[position] != INVALID_CHAR);
        size_t isInvalid = ((*currentChanges)[position] == INVALID_CHAR);
        if (isInvalid) return rootSeq[position];
        return (*currentChanges)[position];

        // ALPHACHAR returnChar =  isChanged ? (*currentChanges)[position] : rootSeq[position];
        // return returnChar;
    }

    

    void handleEvent(const int nodeId, const size_t position, const ALPHACHAR change,
                     const vector<size_t> &rateCategories, const stochasticProcess *sp,
                    sequence &rootSeq) {
        // std::cout << "Change in node: " << nodeId << "\n";
        // std::cout << position << "->" << change << "\n";
        if (_substitutionVec[nodeId] == nullptr) {
            _substitutionVec[nodeId] = std::make_unique<changeMap>(rootSeq.seqLen(), INVALID_CHAR);
        }

        ALPHACHAR previousChar = getCharacter(nodeId,position,rootSeq);

        MDOUBLE previousQii = sp->Qij(previousChar, previousChar);
		MDOUBLE newQii = sp->Qij(change, change);
		updateReactantsSum(-previousQii, sp->rates(rateCategories[position]));
        updateReactantsSum(newQii, sp->rates(rateCategories[position]));


        (*_substitutionVec[nodeId])[position] = change;
        rootSeq[position] = change;
    }

    void dumpSubstitutionLog(const int nodeId) {
        if (nodeId == 0) return;
        if (_substitutionVec[nodeId] == nullptr) return;
        std::stringstream ss;
        std::stringstream changeString;

        ss << "/home/elyawy/temp/" << nodeId << ".sfasta";
        std::ofstream o(ss.str());

        for (auto &changes: *_substitutionVec[nodeId]) {
            // changeString << changes.first << "-" << changes.second << "\n";
            changeString << changes << "\n";
        }
        o << changeString.str();

    }


    void undoSubs(int fromNode, sequence &rootSeq, 
                  const vector<size_t> &rateCategories, const stochasticProcess *sp) {
        if ((_substitutionVec[fromNode] == nullptr)) {
            errorMsg::reportError("Trying to reach removed pointer!");
        }
        auto nodeChangeMap = getChangeMap(fromNode);


        for (size_t currentSite; currentSite < rootSeq.seqLen(); ++currentSite) {
            if ((*nodeChangeMap)[currentSite] == INVALID_CHAR) continue;
            ALPHACHAR currentChar = rootSeq[currentSite];
            ALPHACHAR restoredChar = (*nodeChangeMap)[currentSite];

            MDOUBLE oldFreq = sp->Qij(currentChar, currentChar);
            MDOUBLE newFreq = sp->Qij(restoredChar, restoredChar);
            
            rootSeq[currentSite] = restoredChar;

            updateReactantsSum(-oldFreq, sp->rates(rateCategories[currentSite]));
            updateReactantsSum(newFreq, sp->rates(rateCategories[currentSite]));

	    }


    }



    

    // ALPHACHAR getChange(const int nodeId, const size_t position) {
    //     if ((_substitutionVec[nodeId] == nullptr)) {
    //         errorMsg::reportError("Trying to reach removed pointer!");
    //         return -1;
    //     }

    //     return (*_substitutionVec[nodeId])[position];
    // }

    bool isEmpty(const int nodeId) {

        bool checkEmpty = _substitutionVec[nodeId] == nullptr;
        return checkEmpty;
    }


    std::unique_ptr<changeMap> getChangeMap(const int nodeId) {
        if ((_substitutionVec[nodeId] == nullptr)) {
            std::cout << nodeId << " <- this node is null\n";
            errorMsg::reportError("Trying to reach removed pointer!");
        }

        return std::move(_substitutionVec[nodeId]);
    }

    void printSubManager() {
        std::cout << "printing subs...\n";
        for(auto &node: _substitutionVec) {
            // for(auto &changes: *node){
            //     std::cout << changes.first << "->" << changes.second << ", ";
            // }
            for(auto &changes: *node){
                std::cout << changes << ", ";
            }

            std::cout << "\n";
        }
    }

    void clear() {
        size_t containerSize = _substitutionVec.size();
        _substitutionVec.clear();
        _substitutionVec.resize(containerSize);
        _sumOfReactantsXRates = 0.0;
    }



    ~substitutionManager() {};
};
