#ifndef _MSA
#define _MSA

#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <map>
#include <iostream>
#include <fstream>

#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/sequenceContainer.h"

#include "Sequence.h"


using namespace std;

class MSA
{
public:
	using iteratorType = std::list<SuperSequence::columnContainer>::iterator;

	// MSA(const vector<string> & seqArray) : _originalAlignedSeqs(seqArray), _numberOfSequences(seqArray.size()) {};

    static std::vector<MSA> generateMSAs(const std::vector<BlockMap> &blockmaps, tree::nodeP rootNode) {
        std::vector<MSA> msas;

        for(auto &blockmap: blockmaps) {
            msas.push_back(MSA(blockmap, rootNode));
        }

        return msas;
    }

    MSA (const BlockMap &blockmap,const tree::nodeP rootNode) {
        size_t sequenceSize = std::get<static_cast<int>(BLOCKLIST::LENGTH)>(blockmap.at(rootNode->id()))-1;
        // std::cout << sequenceSize << "\n";

        SuperSequence superSequence(sequenceSize, rootNode->getNumberLeaves());
        Sequence rootSequence(superSequence, rootNode->getNumberOfSons(), false);
        rootSequence.initSequence();
        std::cout << "initialized root\n";
        // std::cin.get();

        std::vector<Sequence> finalSequences;
        buildMsaRecursively(finalSequences, blockmap, *rootNode, superSequence, rootSequence);
        std::cout << "built scaffold MSA\n";

        std::cin.get();
        fillMSA(finalSequences, superSequence);
    }

    void buildMsaRecursively(std::vector<Sequence> &finalSequences,
                             const BlockMap &blockmap,const tree::TreeNode &parrentNode,
                             SuperSequence &superSequence, Sequence &parentSequence) {
        if (parrentNode.isLeaf()) {
            finalSequences.push_back(parentSequence);
            return;
        }
        for (size_t i = 0; i < parrentNode.getNumberOfSons(); i++) {
            tree::TreeNode* childNode = parrentNode.getSon(i);
            Sequence currentSequence(superSequence, childNode->getNumberOfSons(), childNode->isLeaf());
            auto blocks = std::get<static_cast<int>(BLOCKLIST::BLOCKS)>(blockmap.at(childNode->id()));//simulateAlongBranch(sequences.top().size(), currentNode->dis2father(), nodePosition);
            currentSequence.generateSequence(blocks, parentSequence);
            buildMsaRecursively(finalSequences, blockmap, *childNode, superSequence, currentSequence);
        }
        
    }

    void fillMSA(vector<Sequence> &sequences, SuperSequence &superSeq) {
		_numberOfSequences = superSeq.getNumSequences();
		_msaLength = superSeq.getLeafSequenceLength();
		superSeq.setAbsolutePositions();
		_alignedSequence.resize(_numberOfSequences);

        size_t rowInMSA = 0;
        int totalSize = 0;
		int currentPosition = 0;
		int lastPosition = 0;
		int positionDifference = 0;
		int cumulatedDifference = 1;
        
        for(auto &seq: sequences) {
            auto previousSite = *seq.begin();
            lastPosition = previousSite->absolutePosition;
            if (lastPosition > 0) {
                _alignedSequence[rowInMSA].push_back(-lastPosition);
                totalSize += lastPosition;
            }

            for(auto currentSite=seq.begin() + 1; currentSite!=seq.end(); currentSite++) {
				currentPosition = (*(currentSite))->absolutePosition;
                positionDifference = currentPosition - lastPosition - 1;
                
                if (positionDifference == 0) cumulatedDifference++;
                if (positionDifference > 0) {
                    _alignedSequence[rowInMSA].push_back(cumulatedDifference);
                    _alignedSequence[rowInMSA].push_back(-(positionDifference));
                    totalSize += (cumulatedDifference + positionDifference);
                    cumulatedDifference = 1;
                }

                lastPosition = currentPosition;

            }
			if (cumulatedDifference > 0 && (totalSize != _msaLength)) {
                _alignedSequence[rowInMSA].push_back(cumulatedDifference);
                totalSize += cumulatedDifference;
            }
            if (totalSize != _msaLength) _alignedSequence[rowInMSA].push_back(-(_msaLength - totalSize));
			cumulatedDifference = 1;
            lastPosition = 0;
            totalSize = 0;
			rowInMSA++;
        }
    };


    void fillSubstitutions(std::shared_ptr<sequenceContainer> _seqContainer) {
        _substitutions = _seqContainer;
    }

	MSA(size_t numSequences, size_t msaLength): _numberOfSequences(numSequences), _msaLength(msaLength) {};

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

	void printIndels() {

		for (auto const &sequence: _alignedSequence)
		{
            for (auto const &site: sequence) {
                 std::cout << site << " ";     //std::bitset<8>(column);
            }
            std::cout << std::endl;

		}
	}


    std::string generateMsaString() {
        std::stringstream msaString;
        for (size_t row = 0; row < _numberOfSequences; row++) {
            int passedSeq = 0;
            int id = _substitutions->placeToId(row);
            // std::cout << id << "\n";
            msaString << ">" << _substitutions->name(id) << "\n";
            std::string currentSeq = (*_substitutions)[id].toString();
            if (_alignedSequence.empty()) {
                msaString << currentSeq;
                msaString << "\n";
                continue;
            }
            // std::cout << currentSeq << "\n";
            for (size_t col = 0; col < _alignedSequence[row].size(); col++) {
                int strSize = _alignedSequence[row][col];

                if (strSize < 0) {
                    strSize = -strSize;
                    msaString << std::string(strSize, '-');
                } else {
                    msaString << currentSeq.substr(passedSeq, strSize);
                }
                // std::cout << strSize << " ";

                passedSeq += strSize;
            }
            msaString << "\n";
            
        }
        return msaString.str();
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

    std::vector<vector<int>> getMSAVec() {return _alignedSequence;}
	
	~MSA() {
		// delete _originalAlignedSeqs;
	}

private:
	size_t _numberOfSequences; // NUMBER OF SEQUENCES IN THE MSA
    size_t _msaLength; // Length of the MSA
    std::shared_ptr<sequenceContainer> _substitutions;

	SuperSequence* _originalAlignedSeqs; //The aligned sequences

	std::vector<vector<int>> _alignedSequence;

};
#endif
