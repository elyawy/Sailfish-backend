#ifndef _MSA
#define _MSA

#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <map>
#include <iostream>

#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/sequenceContainer.h"

#include "Sequence.h"


using namespace std;

class MSA
{
public:
	using iteratorType = std::list<SuperSequence::columnContainer>::iterator;
    using BlockMap = std::map<std::string, BlockTree>;

	// MSA(const vector<string> & seqArray) : _originalAlignedSeqs(seqArray), _numberOfSequences(seqArray.size()) {};

    static std::vector<MSA> generateMSAs(std::vector<BlockMap> &blockmaps, tree::nodeP rootNode) {
        std::vector<MSA> msas;

        for(auto &blockmap: blockmaps) {
            msas.push_back(MSA(blockmap, rootNode));
        }

        return msas;
    }

    MSA (BlockMap &blockmap, tree::nodeP rootNode) {
        size_t sequenceSize = blockmap[rootNode->name()].length()-1;
        std::stack<tree::nodeP> nodes;
        nodes.push(rootNode);
        tree::nodeP currentNode = nodes.top();

        SuperSequence superSequence(sequenceSize, rootNode->getNumberLeaves());
        Sequence rootSequence(superSequence, currentNode->getNumberOfSons(), false);
        rootSequence.initSequence();

        std::stack<Sequence> sequences;
        sequences.push(rootSequence);
        std::vector<Sequence> finalSequences;

        size_t nodePosition = 0;
        while (!nodes.empty()) {
            nodes.pop();
            if (!currentNode->isLeaf()) {
                for (auto node: currentNode->getSons()) {
                    nodes.push(node);
                }
            } else {
                finalSequences.push_back(sequences.top());
                sequences.pop();
                while (!sequences.empty() && sequences.top().getReferenceCount() == 0) sequences.pop();
                if (nodes.empty()) break;
            }
            currentNode = nodes.top();
            auto blocks = blockmap[currentNode->name()];//simulateAlongBranch(sequences.top().size(), currentNode->dis2father(), nodePosition);
            Sequence currentSequence(superSequence, currentNode->getNumberOfSons(), currentNode->isLeaf());
            currentSequence.generateSequence(blocks, sequences.top());
            // std::cout << "block length " << blocks.length() - 1 << "\n";
            sequences.push(currentSequence);

            ++nodePosition;
        }
        fillMSA(finalSequences, superSequence);
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


    void fillSubstitutions(sequenceContainer* _seqContainer) {
        _substitutions = _seqContainer;
    }


	MSA(){
		_numberOfSequences = 0;
		_msaLength = 0;
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

	void printIndels() {

		for (auto const &sequence: _alignedSequence)
		{
            for (auto const &site: sequence) {
                 std::cout << site << " ";     //std::bitset<8>(column);
            }
            std::cout << std::endl;

		}
	}

    void printFullMsa() {

        for (size_t row = 0; row < _numberOfSequences; row++) {
            int passedSeq = 0;
            int id = _substitutions->placeToId(row);
            // std::cout << id << "\n";
            std::string currentSeq = (*_substitutions)[id].toString();

            // std::cout << currentSeq << "\n";
            for (size_t col = 0; col < _alignedSequence[row].size(); col++) {
                int strSize = _alignedSequence[row][col];

                if (strSize < 0) {
                    strSize = -strSize;
                    std::cout << std::string(strSize, '-');
                } else {
                    std::cout << currentSeq.substr(passedSeq, strSize);
                }
                // std::cout << strSize << " ";

                passedSeq += strSize;
            }
            std::cout << "\n";
            
        }
    
	}



    void writeMSA() {
        
    }

    std::vector<vector<int>> getMSAVec() {return _alignedSequence;}
	
	~MSA() {
		// delete _originalAlignedSeqs;
	}

private:

	SuperSequence* _originalAlignedSeqs; //The aligned sequences

    sequenceContainer* _substitutions;
	std::vector<vector<int>> _alignedSequence;
    size_t _msaLength; // Length of the MSA
	size_t _numberOfSequences; // NUMBER OF SEQUENCES IN THE MSA

};
#endif
