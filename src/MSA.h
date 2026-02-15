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

#include "SimulationContext.h"
#include "Sequence.h"
#include "Event.h"



typedef std::vector<std::vector<int>> SparseMSA;

template<typename RngType = std::mt19937_64>
class MSA
{
public:
	using iteratorType = std::list<SuperSequence::columnContainer>::iterator;


    MSA<RngType>(EventMap &eventmap, SimulationContext<RngType> &simContext): 
        _simContext(simContext),  _numberOfSequences(simContext.getNumberOfNodesToSave()), _msaLength(0) {

        tree::nodeP rootNode = simContext.getRootNode();
        const std::vector<bool>& nodesToSave = simContext.getNodesToSave();
        // Retrieve root sequence size from dummy insertion event in event map
        size_t sequenceSize = eventmap.at(rootNode->id())[0].length;

        // DISPATCH based on protocol setting
        auto* protocol = simContext.getProtocol();
        if (protocol->getIndelRateModel() == IndelRateModel::SIMPLE) {
            buildMsa<BlockTree>(eventmap, rootNode, nodesToSave, sequenceSize);
        } else {
            buildMsa<BlockTreeWithRates>(eventmap, rootNode, nodesToSave, sequenceSize);
        }

    }


    template<typename BlockTreeType>
    void buildMsa(EventMap &eventmap, tree::nodeP rootNode,
                  const std::vector<bool>& nodesToSave, size_t sequenceSize) {
        // Create SuperSequence with the appropriate BlockTree type
        SuperSequence<RngType, BlockTreeType> superSequence(sequenceSize, _numberOfSequences, _simContext);
        
        Sequence rootSequence(superSequence, nodesToSave[rootNode->id()]);
        rootSequence.initSequence();

        std::vector<CompressedSequence> finalSequences;
        finalSequences.reserve(_numberOfSequences);

        buildMsaRecursively(finalSequences, eventmap, *rootNode, 
                           superSequence, rootSequence, nodesToSave);
        eventmap.clear();
        
        fillMSA(finalSequences, superSequence);
    }


    template<typename BlockTreeType>
    void buildMsaRecursively(std::vector<CompressedSequence> &finalSequences,
                             EventMap &eventmap, const tree::TreeNode &parrentNode,
                             SuperSequence<RngType, BlockTreeType> &superSequence, const Sequence& parentSequence, 
                             const std::vector<bool>& nodesToSave) {
        if ((nodesToSave)[parrentNode.id()]) finalSequences.emplace_back(parentSequence.compress());
        if (parrentNode.isLeaf()) return;

        for (size_t i = 0; i < parrentNode.getNumberOfSons(); i++) {
            tree::TreeNode* childNode = parrentNode.getSon(i);
            Sequence currentSequence(superSequence, nodesToSave[childNode->id()]);

            auto events = eventmap.at(childNode->id());
            currentSequence.generateSequence(events, &parentSequence);
            buildMsaRecursively(finalSequences, eventmap, *childNode, superSequence, currentSequence, nodesToSave);
        }
        
    }

    template<typename BlockTreeType>
    void fillMSA(std::vector<CompressedSequence> &sequences, SuperSequence<BlockTreeType> &superSeq) {
		_numberOfSequences = superSeq.getNumSequences();
		_msaLength = superSeq.getMsaSequenceLength();
		superSeq.setAbsolutePositions();
		// create AlignedSequenceMap shared pointer, can be shared with SubstitutionSimulator
        _alignedSequenceMap = std::make_shared<SparseMSA>();
        auto & _alignedSequence = *_alignedSequenceMap;
        _alignedSequence.resize(_numberOfSequences, std::vector<int>());


        int totalSize = 0;
		int currentPosition = 0;
		int lastPosition = 0;
		int positionDifference = 0;
		int cumulatedDifference = 1;
        
        for(size_t i = 0; i < sequences.size(); ++i) {
            Sequence seq(sequences[i], superSeq);
            auto& currentSequence = _alignedSequence[i];
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


    void fillSubstitutions(std::shared_ptr<const SparseSequenceContainer> _seqContainer) {
        _substitutions = _seqContainer;
    }


	MSA<RngType>(size_t msaLength, SimulationContext<RngType> &simContext): 
        _simContext(simContext), _msaLength(msaLength), 
        _numberOfSequences(simContext.getNumberOfNodesToSave()) {
            // create dummy aligned sequence map without gaps.
            _alignedSequenceMap = std::make_shared<SparseMSA>();
            auto & _alignedSequence = *_alignedSequenceMap;
            _alignedSequence.resize(_numberOfSequences, std::vector<int>());
            for (size_t i = 0; i < _numberOfSequences; i++) {
                _alignedSequence[i].push_back(msaLength);
            }
        };

	MSA<RngType>(const MSA<RngType> &msa){
        _simContext = msa._simContext;
		_numberOfSequences = msa._numberOfSequences;
		_msaLength = msa._msaLength;
        _alignedSequenceMap = msa._alignedSequenceMap;
	};

	int getMSAlength() const {return _msaLength;}
	int getNumberOfSequences() const {return _numberOfSequences;} 

	void printMSAInfo() {
		std::cout << _numberOfSequences << "x" << _alignedSequenceMap->size() << "\n";
		std::cout << _msaLength << "\n";
	}





    std::string generateMsaRowString(size_t row) {
        std::string msaString;
        msaString.reserve(_msaLength + 1);

        msaString.append(">");
        msaString.append(_simContext.getNodeToSaveNames()[row]);
        msaString.append("\n");

        std::string currentSeq;

        if (_substitutions == nullptr) {
            currentSeq = "";
        } else {
            currentSeq = (*_substitutions)[row];
        }

        int passedSeq = 0;
        const auto& alignedSeqRow = _alignedSequenceMap->at(row);
        for (size_t col = 0; col < alignedSeqRow.size(); col++) {
            int strSize = alignedSeqRow[col];
            if (strSize < 0) {
                strSize = -strSize;
                msaString.append(strSize, '-');
            } else {
                if (_substitutions == nullptr) {
                    msaString.append(strSize, 'A');
                } else {
                    msaString.append(currentSeq.substr(passedSeq, strSize));
                }
                passedSeq += strSize;
            }
        }
        msaString.append("\n");

        return msaString;
    }

    // print the Msa to console line by line, to avoid creating a large string in memory that represents the whole MSA.
    void printFullMsa() {
        for (size_t row = 0; row < _numberOfSequences; row++) {
            std::cout << generateMsaRowString(row);
        }
    }

    // overload "<<" operator to print MSA to output stream
    friend std::ostream& operator<<(std::ostream& os, MSA<RngType>& msa) {
        for (size_t row = 0; row < msa._numberOfSequences; row++) {
            os << msa.generateMsaRowString(row);
        }
        return os;
    }


    void writeFullMsa(const char * filePath) {
        ofstream msafile (filePath);
        if (msafile.is_open()) {
            for (size_t row = 0; row < _numberOfSequences; row++) {
                msafile << generateMsaRowString(row);
            }
            msafile.close();
        }
        else cout << "Unable to open file";
    }

    std::shared_ptr<const SparseMSA> getSparseMSA() {return _alignedSequenceMap;}

	
	~MSA() {
	}

private:
    const SimulationContext<RngType>& _simContext;
	size_t _numberOfSequences; // NUMBER OF SEQUENCES IN THE FINAL MSA
    size_t _msaLength; // Length of the MSA
    std::shared_ptr<const SparseSequenceContainer> _substitutions;

	std::shared_ptr<SparseMSA> _alignedSequenceMap;
    // IDs of sequences to save are accessed via _simContext->getNodesToSaveIndices()
    // Names of sequences to save are accessed via _simContext->getNodeToSaveNames()


};
#endif
