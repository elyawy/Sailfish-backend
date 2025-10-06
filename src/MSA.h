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


using namespace std;

class MSA
{
public:
	using iteratorType = std::list<SuperSequence::columnContainer>::iterator;

	// MSA(const vector<string> & seqArray) : _originalAlignedSeqs(seqArray), _numberOfSequences(seqArray.size()) {};

    static std::vector<MSA> generateMSAs(const std::vector<BlockMap> &blockmaps, tree::nodeP rootNode,
                                        const std::vector<bool>& nodesToSave) {
        std::vector<MSA> msas;

        for(auto &blockmap: blockmaps) {
            msas.push_back(MSA(blockmap, rootNode, nodesToSave));
        }

        return msas;
    }

    MSA (const BlockMap &blockmap,const tree::nodeP rootNode, const std::vector<bool>& nodesToSave) {
        size_t sequenceSize = std::get<static_cast<int>(BLOCKLIST::LENGTH)>(blockmap.at(rootNode->id()))-1;

        size_t numberOfSeqs = 0;
        _sequencesToSave.clear();
        for (size_t i=0; i < (nodesToSave).size(); i++) {
            numberOfSeqs += (nodesToSave)[i];
            if ((nodesToSave)[i]) _sequencesToSave.push_back(i);
        }
        _numberOfSequences = numberOfSeqs;
        // SuperSequence superSequence(sequenceSize, rootNode->getNumberLeaves());
        SuperSequence superSequence(sequenceSize, numberOfSeqs);
        std::shared_ptr<Sequence> rootSequence =
                        std::make_shared<Sequence>(superSequence, nodesToSave[rootNode->id()], rootNode->id());
        rootSequence->initSequence();
        // std::cin.get();
        std::vector<std::shared_ptr<Sequence>> finalSequences;
        // std::vector<Sequence> finalSequences;

        buildMsaRecursively(finalSequences, blockmap, *rootNode, superSequence, rootSequence, nodesToSave);
        fillMSA(finalSequences, superSequence);
    }

    void buildMsaRecursively(std::vector<std::shared_ptr<Sequence>> &finalSequences,
                             const BlockMap &blockmap,const tree::TreeNode &parrentNode,
                             SuperSequence &superSequence, std::shared_ptr<Sequence> parentSequence, 
                             const std::vector<bool>& nodesToSave) {
        if ((nodesToSave)[parrentNode.id()]) finalSequences.push_back(parentSequence);
        if (parrentNode.isLeaf()) return;

        for (size_t i = 0; i < parrentNode.getNumberOfSons(); i++) {
            tree::TreeNode* childNode = parrentNode.getSon(i);
            std::shared_ptr<Sequence> currentSequence = 
                        std::make_shared<Sequence>(superSequence, nodesToSave[childNode->id()], childNode->id());

            auto blocks = std::get<static_cast<int>(BLOCKLIST::BLOCKS)>(blockmap.at(childNode->id()));//simulateAlongBranch(sequences.top().size(), currentNode->dis2father(), nodePosition);
            currentSequence->generateSequence(blocks, *parentSequence);
            buildMsaRecursively(finalSequences, blockmap, *childNode, superSequence, currentSequence, nodesToSave);

        }
        
    }

    void fillMSA(vector<std::shared_ptr<Sequence>> &sequences, SuperSequence &superSeq) {
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
        
        for(auto &seq: sequences) {
            size_t sequenceNodeID = seq->getSequenceNodeID();
            auto& currentSequence = _alignedSequence[sequenceNodeID];
            // if the sequence is only made up of gaps:
            if (seq->size() == 0) {
                currentSequence.push_back(-_msaLength);
                continue;
            }
            auto previousSite = *seq->begin();
            
            lastPosition = previousSite->absolutePosition;
            if (lastPosition > 0) {
                currentSequence.push_back(-lastPosition);
                totalSize += lastPosition;
            }

            for(auto currentSite=seq->begin() + 1; currentSite!=seq->end(); currentSite++) {
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

    void setSubstitutionsFolder(const std::string& substitutionsDir) {
        _substitutionsDir = substitutionsDir;
// std::stoi(entry.path().stem())

        for (const auto& entry : std::filesystem::directory_iterator(_substitutionsDir)) {
            
            _substitutionPaths.push_back(entry);
        }
        // std::cout << _substitutionPaths.size() << "\n";
    }


	MSA(size_t numSequences, size_t msaLength,const std::vector<bool>& nodesToSave): 
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

	void printIndels() {

		for (auto const &sequence: _alignedSequence)
		{
            for (auto const &site: sequence.second) {
                 std::cout << site << " ";     //std::bitset<8>(column);
            }
            std::cout << std::endl;

		}
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


    void writeMsaFromDir(const char * filePath) {
        std::ofstream msafile (filePath);

        if (msafile.is_open()) {
            for (size_t row = 0; row < _numberOfSequences; row++) {
                int passedSeq = 0;
                auto & seqPath = _substitutionPaths[row].path();
                int id = std::stoi(seqPath.stem());

                std::ifstream seqFile(seqPath);
                char currentChar;
                while (seqFile.get(currentChar)) {
                    msafile << currentChar;
                    if (currentChar == '\n') break;;
                }
                
                if (_alignedSequence.empty()) {
                    while (seqFile.get(currentChar)) {
                        msafile << currentChar;
                    }
                    continue;
                }

                for (size_t col = 0; col < _alignedSequence[id].size(); col++) {
                    int strSize = _alignedSequence[id][col];
                    if (strSize < 0) {
                        strSize = -strSize;
                        size_t readCounter = strSize;
                        while (readCounter > 0) {
                            msafile << '-';
                            seqFile.get(currentChar); 
                            readCounter--;
                        }
                    } else {
                        size_t readCounter = strSize;
                        while (readCounter > 0) {
                            seqFile.get(currentChar);
                            msafile << currentChar;
                            readCounter--;
                        }
                    }
                    passedSeq += strSize;
                }
                msafile << '\n'; 
                seqFile.close(); 
                std::filesystem::remove(seqPath); // delete taxa sequence file              
            }
            msafile.close();
        }
        else cout << "Unable to open file";
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
	
	~MSA() {
		// delete _originalAlignedSeqs;
	}

private:
	size_t _numberOfSequences; // NUMBER OF SEQUENCES IN THE MSA
    size_t _msaLength; // Length of the MSA
    std::shared_ptr<sequenceContainer> _substitutions;
    std::string _substitutionsDir;
    std::vector<std::filesystem::directory_entry> _substitutionPaths;

	SuperSequence* _originalAlignedSeqs; //The aligned sequences

	std::unordered_map<size_t, std::vector<int>> _alignedSequence;
    std::vector<size_t> _sequencesToSave;


};
#endif
