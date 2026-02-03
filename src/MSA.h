#ifndef _MSA
#define _MSA

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

#include "Sequence.h"
#include "Event.h"

using namespace std;

class MSA
{
public:
	using iteratorType = std::list<SuperSequence::columnContainer>::iterator;


    MSA (EventMap &eventmap, size_t sequenceSize, const tree::nodeP rootNode, const std::vector<bool>& nodesToSave) {

        size_t numberOfSeqs = 0;
        _sequencesToSave.clear();
        for (size_t i=0; i < (nodesToSave).size(); i++) {
            numberOfSeqs += (nodesToSave)[i];
            if ((nodesToSave)[i]) _sequencesToSave.push_back(i);
        }
        _numberOfSequences = numberOfSeqs;
        // SuperSequence superSequence(sequenceSize, rootNode->getNumberLeaves());
        SuperSequence superSequence(sequenceSize, numberOfSeqs);
        Sequence rootSequence(superSequence, nodesToSave[rootNode->id()], rootNode->id());
        rootSequence.initSequence();

        // std::cin.get();
        std::vector<CompressedSequence> finalSequences;
        finalSequences.reserve(_numberOfSequences);
        // std::vector<Sequence> finalSequences;

        buildMsaRecursively(finalSequences, eventmap, *rootNode, superSequence, rootSequence, nodesToSave);
        eventmap.clear();
        
        fillMSA(finalSequences, superSequence);
    }

    void buildMsaRecursively(std::vector<CompressedSequence> &finalSequences,
                             EventMap &eventmap, const tree::TreeNode &parrentNode,
                             SuperSequence &superSequence, const Sequence& parentSequence, 
                             const std::vector<bool>& nodesToSave) {
        if ((nodesToSave)[parrentNode.id()]) finalSequences.emplace_back(parentSequence.compress());
        if (parrentNode.isLeaf()) return;

        for (size_t i = 0; i < parrentNode.getNumberOfSons(); i++) {
            tree::TreeNode* childNode = parrentNode.getSon(i);
            Sequence currentSequence(superSequence, nodesToSave[childNode->id()], childNode->id());

            auto events = eventmap.at(childNode->id());
            currentSequence.generateSequence(events, &parentSequence);
            buildMsaRecursively(finalSequences, eventmap, *childNode, superSequence, currentSequence, nodesToSave);
        }
        
    }

    void fillMSA(std::vector<CompressedSequence> &sequences, SuperSequence &superSeq) {
		_numberOfSequences = superSeq.getNumSequences();
		_msaLength = superSeq.getMsaSequenceLength();
		superSeq.setAbsolutePositions();
		// _alignedSequence.resize(_numberOfSequences);
        _alignedSequence.reserve(_numberOfSequences);

        // size_t rowInMSA = 0;
        int totalSize = 0;
		int currentPosition = 0;
		int lastPosition = 0;
		int positionDifference = 0;
		int cumulatedDifference = 1;
        
        for(auto &_seq: sequences) {
            Sequence seq(_seq,superSeq);
            size_t sequenceNodeID = seq.getSequenceNodeID();
            auto& currentSequence = _alignedSequence[sequenceNodeID];
            // if the sequence is only made up of gaps:
            if (seq.size() == 0) {
                currentSequence.push_back(-_msaLength);
                continue;
            }
            auto previousSite = *seq.begin();
            
            lastPosition = previousSite->absolutePosition;
            if (lastPosition > 0) {
                currentSequence.push_back(-lastPosition);
                totalSize += lastPosition;
            }

            for(auto currentSite=seq.begin() + 1; currentSite!=seq.end(); currentSite++) {
				currentPosition = (*(currentSite))->absolutePosition;
                positionDifference = currentPosition - lastPosition - 1;

                if (positionDifference == 0) cumulatedDifference++;
                if (positionDifference > 0) {
                    currentSequence.push_back(cumulatedDifference);
                    currentSequence.push_back(-(positionDifference));
                    totalSize += (cumulatedDifference + positionDifference);
                    cumulatedDifference = 1;
                }
                if (totalSize > _msaLength) errorMsg::reportError("sequence lengths mismatch in fillMSA");


                lastPosition = currentPosition;

            }
			if (cumulatedDifference > 0 && (totalSize != _msaLength)) {
                currentSequence.push_back(cumulatedDifference);
                totalSize += cumulatedDifference;
            }
            if (totalSize < _msaLength) currentSequence.push_back(-(_msaLength - totalSize));
			cumulatedDifference = 1;
            lastPosition = 0;
            totalSize = 0;
			// rowInMSA++;
        }
    };

    static MSA msaFromSequences(vector<Sequence> &sequences, SuperSequence &superSeq) {
        std::vector<bool> seqsToSave;
        for (size_t i = 0; i < sequences.size(); i++)
        {
            seqsToSave.push_back(true);
        }
        

        MSA msa(sequences.size(), 1, seqsToSave);
		msa._numberOfSequences = superSeq.getNumSequences();
		msa._msaLength = superSeq.getMsaSequenceLength();
		superSeq.setAbsolutePositions();
		// _alignedSequence.resize(_numberOfSequences);
        msa._alignedSequence.reserve(msa._numberOfSequences);
        
        // size_t rowInMSA = 0;
        int totalSize = 0;
		int currentPosition = 0;
		int lastPosition = 0;
		int positionDifference = 0;
		int cumulatedDifference = 1;
        
        for(auto &seq: sequences) {
            size_t sequenceNodeID = seq.getSequenceNodeID();
            // if the sequence is only made up of gaps:
            if (seq.size() == 0) {
                msa._alignedSequence[sequenceNodeID].push_back(-msa._msaLength);
                continue;
            }
            // seq.printSequence();

            auto previousSite = *seq.begin();
            
            lastPosition = previousSite->absolutePosition;
            if (lastPosition > 0) {
                msa._alignedSequence[sequenceNodeID].push_back(-lastPosition);
                totalSize += lastPosition;
            }

            for(auto currentSite=seq.begin() + 1; currentSite!=seq.end(); currentSite++) {
				currentPosition = (*(currentSite))->absolutePosition;
                positionDifference = currentPosition - lastPosition - 1;

                if (positionDifference == 0) cumulatedDifference++;
                if (positionDifference > 0) {
                    msa._alignedSequence[sequenceNodeID].push_back(cumulatedDifference);
                    msa._alignedSequence[sequenceNodeID].push_back(-(positionDifference));
                    totalSize += (cumulatedDifference + positionDifference);
                    cumulatedDifference = 1;
                }
                if (totalSize > msa._msaLength) errorMsg::reportError("sequence lengths mismatch in fillMSA");


                lastPosition = currentPosition;

            }
			if (cumulatedDifference > 0 && (totalSize != msa._msaLength)) {
                msa._alignedSequence[sequenceNodeID].push_back(cumulatedDifference);
                totalSize += cumulatedDifference;
            }
            if (totalSize < msa._msaLength) msa._alignedSequence[sequenceNodeID].push_back(-(msa._msaLength - totalSize));
			cumulatedDifference = 1;
            lastPosition = 0;
            totalSize = 0;
			// rowInMSA++;
        }
        return msa;
    };


    void fillSubstitutions(std::shared_ptr<sequenceContainer> _seqContainer) {
        _substitutions = _seqContainer;
    }


	MSA(size_t numSequences, size_t msaLength, const std::vector<bool>& nodesToSave): 
        _numberOfSequences(numSequences), _msaLength(msaLength) {
        _sequencesToSave.clear();
        for (size_t i=0; i < (nodesToSave).size(); i++) {
            if ((nodesToSave)[i]) _sequencesToSave.push_back(i);
        }
    };

	MSA(const MSA &msa){
		_numberOfSequences = msa._numberOfSequences;
		_msaLength = msa._msaLength;
		_alignedSequence = msa._alignedSequence;
	};

	int getMSAlength() const {return _msaLength;}
	int getNumberOfSequences() const {return _numberOfSequences;} 

	void printMSAInfo() {
		std::cout << _numberOfSequences << "x" << _alignedSequence.size() << "\n";
		std::cout << _msaLength << "\n";
	}


    std::string generateMsaStringWithoutSubs() {
        // std::stringstream msaString;
        std::string msaString;
        msaString.reserve((_msaLength+1)*_numberOfSequences);

        for (auto id: _sequencesToSave) {
            const auto& alignedSeqRow = _alignedSequence[id];  // Single hash lookup

            for (size_t col = 0; col < alignedSeqRow.size(); col++) {
                int strSize = alignedSeqRow[col];
                if (strSize < 0) {
                    strSize = -strSize;
                    msaString.append(strSize, '-');
                } else {
                    msaString.append(strSize, 'A');
                }
            }
            msaString.append("\n");
        }
        return msaString;
    }


    std::string generateMsaString() {
        if (_substitutions == nullptr) return generateMsaStringWithoutSubs();
        std::string msaString;
        msaString.reserve((_msaLength+1)*_numberOfSequences);
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
            const auto& alignedSeqRow = _alignedSequence[id];  // Single hash lookup
            for (size_t col = 0; col < alignedSeqRow.size(); col++) {
                int strSize = alignedSeqRow[col];
                if (strSize < 0) {
                    strSize = -strSize;
                    msaString.append(strSize, '-');

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

    std::unordered_map<size_t, std::vector<int>> getMSAVec() {return _alignedSequence;}

    const std::unordered_map<size_t, std::vector<int>>& getAlignedSequence() const {
        return _alignedSequence;
    }
	
	~MSA() {
	}

private:
	size_t _numberOfSequences; // NUMBER OF SEQUENCES IN THE MSA
    size_t _msaLength; // Length of the MSA
    std::shared_ptr<sequenceContainer> _substitutions;

	std::unordered_map<size_t, std::vector<int>> _alignedSequence;
    std::vector<size_t> _sequencesToSave;


};
#endif
