#ifndef _MSA_FIXED
#define _MSA_FIXED

#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/sequenceContainer.h"

#include "IteratorSequence.h"
#include "FixedList.h"

using namespace std;

class MsaFixed
{
public:
    using iteratorType = FixedList::iterator;


    MsaFixed(const BlockMap &blockmap, const tree::nodeP rootNode, 
         const std::vector<bool>& nodesToSave) {
        size_t sequenceSize = std::get<static_cast<int>(BLOCKLIST::LENGTH)>(blockmap.at(rootNode->id()))-1;

        size_t numberOfSeqs = 0;
        _sequencesToSave.clear();
        for (size_t i=0; i < (nodesToSave).size(); i++) {
            numberOfSeqs += (nodesToSave)[i];
            if ((nodesToSave)[i]) _sequencesToSave.push_back(i);
        }
        _numberOfSequences = numberOfSeqs;

        // Calculate exact size needed by FixedList
        size_t totalInsertions = 0;
        for (const auto& [nodeId, blockTuple] : blockmap) {
            const BlockList& blocks = std::get<static_cast<int>(BLOCKLIST::BLOCKS)>(blockTuple);
            for (const auto& block : blocks) {
                totalInsertions += block[static_cast<int>(BLOCK::INSERTION)];
            }
        }

        size_t exactSize = sequenceSize + totalInsertions + 1;

        FixedList fixedList(exactSize);
        fixedList.initialize(sequenceSize);

        IteratorSequence rootSequence(fixedList, nodesToSave[rootNode->id()], rootNode->id());
        rootSequence.initSequence();

        std::vector<IteratorSequence> finalSequences;
        buildMsaRecursivelyFixed(finalSequences, blockmap, *rootNode, fixedList, 
                                rootSequence, nodesToSave);
        fillMSAFixed(finalSequences, fixedList);
    }

    void buildMsaRecursivelyFixed(std::vector<IteratorSequence> &finalSequences,
                                   const BlockMap &blockmap,
                                   const tree::TreeNode &parentNode,
                                   FixedList &fixedList,
                                   IteratorSequence &parentSequence, 
                                   const std::vector<bool>& nodesToSave) {
        if ((nodesToSave)[parentNode.id()]) finalSequences.push_back(parentSequence);
        if (parentNode.isLeaf()) return;

        for (size_t i = 0; i < parentNode.getNumberOfSons(); i++) {
            tree::TreeNode* childNode = parentNode.getSon(i);
            IteratorSequence currentSequence(fixedList, nodesToSave[childNode->id()], childNode->id());
            auto blocks = std::get<static_cast<int>(BLOCKLIST::BLOCKS)>(blockmap.at(childNode->id()));
            currentSequence.generateSequence(blocks, parentSequence);
            buildMsaRecursivelyFixed(finalSequences, blockmap, *childNode, fixedList, 
                                      currentSequence, nodesToSave);
        }
    }

    void fillMSAFixed(vector<IteratorSequence> &sequences, FixedList &fixedList) {
        fixedList.setAbsolutePositions();
        _msaLength = fixedList.getMsaSequenceLength();
        std::cout << "MSaFixedLength=" << _msaLength << "\n";
        std::cout << "FixedList=";
        fixedList.printSequence();

        _alignedSequence.reserve(_numberOfSequences);
        
        int totalSize = 0;
        int currentPosition = 0;
        int lastPosition = 0;
        int positionDifference = 0;
        int cumulatedDifference = 1;
        
        for(auto &seq: sequences) {
            size_t sequenceNodeID = seq.getSequenceNodeID();
            std::cout << "Sequence=" << sequenceNodeID  << ":\n";
            seq.printSequence();
            // if the sequence is only made up of gaps:
            if (seq.size() == 0) {
                _alignedSequence[sequenceNodeID].push_back(-_msaLength);
                continue;
            }

            auto previousSite = *(seq.begin());

            std::cout << *previousSite << "\n";
            lastPosition = fixedList.getAbsolutePosition(previousSite);
            if (lastPosition > 0) {
                _alignedSequence[sequenceNodeID].push_back(-lastPosition);
                totalSize += lastPosition;
            }
            std::cout << "currentPositions=" ;

            for(auto currentSite=seq.begin() + 1; currentSite!=seq.end(); currentSite++) {

                currentPosition = fixedList.getAbsolutePosition(*currentSite);
                positionDifference = currentPosition - lastPosition - 1;
                std::cout << currentPosition << " ";

                if (positionDifference == 0) cumulatedDifference++;
                if (positionDifference > 0) {
                    _alignedSequence[sequenceNodeID].push_back(cumulatedDifference);
                    _alignedSequence[sequenceNodeID].push_back(-(positionDifference));
                    totalSize += (cumulatedDifference + positionDifference);
                    cumulatedDifference = 1;
                }
                if (totalSize > _msaLength) errorMsg::reportError("sequence lengths mismatch in fillMSAFixed");

                lastPosition = currentPosition;
            }
            
            if (cumulatedDifference > 0 && (totalSize != _msaLength)) {
                _alignedSequence[sequenceNodeID].push_back(cumulatedDifference);
                totalSize += cumulatedDifference;
            }
            if (totalSize < _msaLength) _alignedSequence[sequenceNodeID].push_back(-(_msaLength - totalSize));
            cumulatedDifference = 1;
            lastPosition = 0;
            totalSize = 0;
        }
    }

    void fillSubstitutions(std::shared_ptr<sequenceContainer> _seqContainer) {
        _substitutions = _seqContainer;
    }

    MsaFixed(size_t numSequences, size_t msaLength, const std::vector<bool>& nodesToSave): 
        _numberOfSequences(numSequences), _msaLength(msaLength) {
        _sequencesToSave.clear();
        for (size_t i=0; i < (nodesToSave).size(); i++) {
            if ((nodesToSave)[i]) _sequencesToSave.push_back(i);
        }
    }

    MsaFixed(const MsaFixed &msa) {
        _numberOfSequences = msa._numberOfSequences;
        _msaLength = msa._msaLength;
        _alignedSequence = msa._alignedSequence;
    }

    int getMSAlength() const { return _msaLength; }
    int getNumberOfSequences() const { return _numberOfSequences; }

    void printMSAInfo() {
        std::cout << _numberOfSequences << "x" << _alignedSequence.size() << "\n";
        std::cout << _msaLength << "\n";
    }

    void printIndels() {
        for (auto const &sequence: _alignedSequence) {
            for (auto const &site: sequence.second) {
                std::cout << site << " ";
            }
            std::cout << std::endl;
        }
    }

    std::string generateMsaStringWithoutSubs() {
        std::stringstream msaString;
        for (auto id: _sequencesToSave) {
            for (size_t col = 0; col < _alignedSequence[id].size(); col++) {
                int strSize = _alignedSequence[id][col];
                if (strSize < 0) {
                    strSize = -strSize;
                    msaString << std::string(strSize, '-');
                } else {
                    msaString << std::string(strSize, 'A');
                }
            }
            msaString << "\n";
        }
        return msaString.str();
    }

    std::string generateMsaString() {
        if (_substitutions == nullptr) return generateMsaStringWithoutSubs();
        std::string msaString;
        msaString.reserve((_msaLength+256)*_numberOfSequences);
        for (size_t row = 0; row < _numberOfSequences; row++) {
            int passedSeq = 0;
            int id = _substitutions->placeToId(row);
            msaString.append(">");
            msaString.append(_substitutions->name(id));
            msaString.append("\n");
            std::string currentSeq = (*_substitutions)[id].toString();
            if (_alignedSequence.empty()) {
                msaString.append(currentSeq);
                msaString.append("\n");
                continue;
            }
            for (size_t col = 0; col < _alignedSequence[id].size(); col++) {
                int strSize = _alignedSequence[id][col];
                if (strSize < 0) {
                    strSize = -strSize;
                    msaString.append(std::string(strSize, '-'));
                } else {
                    msaString.append(currentSeq.substr(passedSeq, strSize));
                }
                passedSeq += strSize;
            }
            msaString.append("\n");
        }
        return msaString;
    }

    void printFullMsa() {
        std::cout << generateMsaString();
    }

    void writeFullMsa(const char * filePath) {
        ofstream msafile (filePath);
        if (msafile.is_open()) {
            msafile << generateMsaString();
            msafile.close();
        }
        else cout << "Unable to open file";
    }

    std::unordered_map<size_t, std::vector<int>> getMSAVec() { return _alignedSequence; }
    
    ~MsaFixed() {}

private:
    size_t _numberOfSequences;
    size_t _msaLength;
    std::shared_ptr<sequenceContainer> _substitutions;
    std::unordered_map<size_t, std::vector<int>> _alignedSequence;
    std::vector<size_t> _sequencesToSave;
};

#endif