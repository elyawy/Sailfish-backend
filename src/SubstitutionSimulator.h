#ifndef ___SUBSTITUTION_SIMULATOR_H
#define ___SUBSTITUTION_SIMULATOR_H

#include "../libs/Phylolib/includes/definitions.h"
#include "../libs/Phylolib/includes/stochasticProcess.h"

#include "SimulationContext.h"
#include "MSA.h"
#include "modelFactory.h"
#include "CategorySampler.h"
#include "BranchTransitionProbabilities.h"


template<typename RngType = std::mt19937_64,size_t AlphabetSize = 4>
class SubstitutionSimulator {
public:
	explicit SubstitutionSimulator(modelFactory& mFac, SimulationContext<RngType>& simContext): 
		_tree(simContext.getTree()), 
		_stochasticProcess(mFac.getStochasticProcess()), _alphabet(mFac.getAlphabet()), 
		_nodesToSave(simContext.getNodesToSave()), _idToRowInMSA(simContext.getIdToSaveIndices()),
		_saveRates(false),
		_rateCategorySampler(mFac.getEffectiveTransitionMatrix(), mFac.getStationaryProbs()),
		_rng(simContext.getRng()),
		_finalMsaPath("") {

		std::vector<MDOUBLE> frequencies;
		for (int j = 0; j < AlphabetSize; ++j) {
			frequencies.push_back(_stochasticProcess->freq(j));
			_charLookup[j] = _alphabet->fromInt(j);
		}

		_frequencySampler = std::make_unique<DiscreteDistribution>(frequencies);
		_simulatedSequences = std::make_unique<SparseSequenceContainer>();
	}

	virtual ~SubstitutionSimulator() {
		if (_outputFile.is_open()) {
			_outputFile.close();
		}
	}

	void initSubstitionSim(modelFactory& mFac) {
		_stochasticProcess = mFac.getStochasticProcess();
		_alphabet = mFac.getAlphabet();

		_rateCategorySampler = CategorySampler(mFac.getEffectiveTransitionMatrix(), mFac.getStationaryProbs());
		_frequencySampler = std::make_unique<DiscreteDistribution>(mFac.getStationaryProbs());
		_simulatedSequences = std::make_unique<SparseSequenceContainer>();

    }


	void setSaveRates(bool saveRates) {
		_saveRates = saveRates;
	}

	void setPerSiteRateCategories(std::shared_ptr<const std::vector<size_t>> rateCategories) {
		_rateCategories = rateCategories;
	}


	void clearRatesVec() { 
		_siteRates.clear();
	}


	std::unique_ptr<SparseSequenceContainer> getSequenceContainer() {
		auto outputSequences = std::move(_simulatedSequences);
		_simulatedSequences = std::make_unique<SparseSequenceContainer>();
		return outputSequences;
	}


	std::vector<double> getSiteRates() { 
		return _siteRates;
	}

	void generateSubstitutionsAlongTree(int seqLength) {

		if (_rateCategories == nullptr) {
			auto newCategories = std::make_shared<std::vector<size_t>>(seqLength);
			for (int h = 0; h < seqLength; h++) {
            	(*newCategories)[h] = _rateCategorySampler.drawSample(_rng);
			}
			_rateCategories = newCategories;
        } else {
			// Use existing categories - just validate size
			if (_rateCategories->size() != seqLength) {
				errorMsg::reportError("Rate categories size mismatch");
			}
		}

		_siteRates.clear();
		if (_saveRates){
			_siteRates.resize(seqLength);
			for (int i = 0; i < seqLength; i++) {
				_siteRates[i] = _stochasticProcess->rates((*_rateCategories)[i]);
			}
		}
		
		sequence rootSequence = generateRootSeq(seqLength);

		if (_nodesToSave[_tree->getRoot()->id()]){ 
			saveSequence(rootSequence);
		}

		mutateSeqRecuresively(rootSequence, _tree->getRoot());
	}

	void mutateSeqRecuresively(const sequence& currentSequence, tree::nodeP currentNode) {
		if (currentNode->isLeaf()) return;
		for (auto &node: currentNode->getSons()) {
			sequence childSeq(currentSequence);
			childSeq.setID(node->id());
			childSeq.setName(node->name());
			mutateSeqAlongBranch(childSeq, node->dis2father());
			if (_nodesToSave[node->id()]) saveSequence(childSeq);
			mutateSeqRecuresively(childSeq, node);
		}
	}

	void setWriteFolder(const std::string &filePath) {
		_finalMsaPath = filePath;
		if (!filePath.empty()) {
			_outputFile.open(filePath, std::ios::out);  // open once
			if (!_outputFile.is_open()) {
				errorMsg::reportError("Could not open file " + filePath + " for writing MSA.");
			}
		}
	}


    void simulateAndWriteSubstitutions(size_t sequenceLength, const std::string& filePath) {
        setWriteFolder(filePath);
        generateSubstitutionsAlongTree(sequenceLength);
    }

    std::shared_ptr<SparseSequenceContainer> simulateSubstitutions(size_t sequenceLength) {
        generateSubstitutionsAlongTree(sequenceLength);

        return getSequenceContainer();
    }


    void setAlignedSequenceMap(MSA<RngType>& msa) {
        _alignedSequenceMap = msa.getSparseMSA();
    }


private:

	sequence generateRootSeq(int seqLength) {
		sequence rootSeq(_alphabet);

		rootSeq.resize(seqLength);

		size_t rootID = _tree->getRoot()->id();
		for (int i = 0; i < seqLength; i++) {
			ALPHACHAR newChar = _frequencySampler->drawSample(_rng) - 1;
			rootSeq[i] = newChar;
		}
		rootSeq.setName(_tree->getRoot()->name());
		rootSeq.setID(_tree->getRoot()->id());
		return rootSeq;

	}

	void mutateSeqAlongBranch(sequence& currentSequence, const MDOUBLE& distToFather) {
		mutateEntireSeq(currentSequence, distToFather);
	}

	void mutateEntireSeq(sequence& currentSequence, const MDOUBLE& branchLength) {
		const int nodeId = currentSequence.id();
		auto& rateCategories = (*_rateCategories);
		BranchTransitionProbabilities<AlphabetSize> cachedPijt(branchLength, *_stochasticProcess);
		size_t actualRowInMSA = _idToRowInMSA[nodeId];
		// Check if this is a leaf we're saving (low memory mode)
		if (_alignedSequenceMap != nullptr && _nodesToSave[nodeId]) {
			const std::vector<int>& gapStructure = _alignedSequenceMap->at(actualRowInMSA);
			
			size_t site = 0;
			_lengthOfCurrentSequence = 0;
			for (int blockSize : gapStructure) {
				if (blockSize < 0) {
					// Gap block - skip these sites
					site += (-blockSize);
					continue;
				}
				// Non-gap block - mutate these sites
				for (int i = 0; i < blockSize; ++i, ++site) {
					ALPHACHAR parentChar = currentSequence[site];
					auto &Pijt = cachedPijt.getDistribution(rateCategories[site], parentChar);
					ALPHACHAR nextChar = Pijt.drawSample(_rng) - 1;
					currentSequence[site] = nextChar;
				}
				_lengthOfCurrentSequence += (blockSize);
			}
		} else {
			// Normal mode - mutate all sites
			for (size_t site = 0; site < currentSequence.seqLen(); ++site) {
				ALPHACHAR parentChar = currentSequence[site];
				auto &Pijt = cachedPijt.getDistribution(rateCategories[site], parentChar);
				ALPHACHAR nextChar = Pijt.drawSample(_rng) - 1;
				currentSequence[site] = nextChar;
			}
		}
	}


	void saveSequence(const sequence &currentSequence) {
		if (_finalMsaPath.size() > 0) {
			saveSequenceToDisk(currentSequence);
			return;
		}
		SparseSequence sparseSeq;
		sparseSeq.reserve(_lengthOfCurrentSequence);
		// populate the sparse sequence with only non-gap characters
		// using the aligned sequence map if available to determine where the gaps are
		if (_alignedSequenceMap != nullptr) {
			size_t actualRowInMSA = _idToRowInMSA[currentSequence.id()];
			const std::vector<int>& gapStructure = _alignedSequenceMap->at(actualRowInMSA);
			size_t site = 0;
			for (int blockSize : gapStructure) {
				if (blockSize < 0) {
					// Gap block - skip these sites
					site += (-blockSize);
				} else {
					// Non-gap block - add these characters to the sparse sequence
					for (int i = 0; i < blockSize; ++i, ++site) {
						// build the string from the char from lookup
						sparseSeq += _charLookup[currentSequence[site]];
					}
				}
			}
		} else {
			for (size_t site = 0; site < currentSequence.seqLen(); ++site) {
				sparseSeq += (_charLookup[currentSequence[site]]);
			}
		}
		_simulatedSequences->push_back(std::move(sparseSeq));
	}

	void saveSequenceToDisk(const sequence &currentSequence) {
		const int nodeId = currentSequence.id();
		size_t actualRowInMSA = _idToRowInMSA[nodeId];
		
		_outputFile << ">" << currentSequence.name() << "\n";
		
		// Get gap structure for this sequence
		if (_alignedSequenceMap != nullptr) {
			const std::vector<int>& gapStructure = _alignedSequenceMap->at(actualRowInMSA);
			size_t site = 0;
			for (int blockSize : gapStructure) {
				if (blockSize < 0) {
					// Gap block - write gaps
					_outputFile << std::string(-blockSize, '-');
				} else {
					// Non-gap block - write sequence characters
					for (int i = 0; i < blockSize; ++i) {
						_outputFile << _charLookup[currentSequence[site++]];
					}
				}
			}
		} else {
			for (size_t site = 0; site < currentSequence.seqLen(); ++site) {
				_outputFile << _charLookup[currentSequence[site]];
			}
		}
		_outputFile << "\n";
	}



	tree* _tree;
	std::shared_ptr<const stochasticProcess> _stochasticProcess;
	const alphabet* _alphabet;

	const std::vector<bool>& _nodesToSave;
	const std::vector<size_t>& _idToRowInMSA;
	bool _saveRates;

	std::shared_ptr<const std::vector<size_t>> _rateCategories = nullptr;
	std::vector<double> _siteRates;
	std::unique_ptr<SparseSequenceContainer> _simulatedSequences;
	std::unique_ptr<DiscreteDistribution> _frequencySampler;

	CategorySampler _rateCategorySampler;
	std::string _finalMsaPath;

	std::array<std::string, AlphabetSize> _charLookup;

	std::shared_ptr<const SparseMSA> _alignedSequenceMap = nullptr;

	RngType &_rng;
	std::ofstream _outputFile;

	size_t _lengthOfCurrentSequence = 0;

};

#endif