#ifndef ___SUBSTITUTION_SIMULATOR_H
#define ___SUBSTITUTION_SIMULATOR_H

#include "../libs/Phylolib/includes/definitions.h"
#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/stochasticProcess.h"
#include "../libs/Phylolib/includes/sequenceContainer.h"

#include "MSA.h"
#include "modelFactory.h"
#include "CategorySampler.h"
#include "BranchTransitionProbabilities.h"

template<typename RngType = std::mt19937_64,size_t AlphabetSize = 4>
class SubstitutionSimulator {
public:
	explicit SubstitutionSimulator(modelFactory& mFac, std::shared_ptr<std::vector<bool>> nodesToSave) : 
		_et(mFac.getTree()), _sp(mFac.getStochasticProcess()), _alph(mFac.getAlphabet()), 
		_nodesToSave(nodesToSave), _saveRates(false),
		_rateCategorySampler(mFac.getEffectiveTransitionMatrix(), mFac.getStationaryProbs()),
		_finalMsaPath("") {
		
		
		std::vector<MDOUBLE> frequencies;
		for (int j = 0; j < AlphabetSize; ++j) {
			frequencies.push_back(_sp->freq(j));
			_charLookup[j] = _alph->fromInt(j);
		}

		_frequencySampler = std::make_unique<DiscreteDistribution>(frequencies);
		_simulatedSequences = std::make_unique<sequenceContainer>();
	}

	virtual ~SubstitutionSimulator() {
		if (_outputFile.is_open()) {
			_outputFile.close();
		}
	}

	//init the substitution simulator with the given model factory
	void initSubstitionSim(modelFactory& mFac) {
		_et = mFac.getTree();
		_sp = mFac.getStochasticProcess();
		_alph = mFac.getAlphabet();

		_rateCategorySampler = CategorySampler(mFac.getEffectiveTransitionMatrix(), mFac.getStationaryProbs());
		_frequencySampler = std::make_unique<DiscreteDistribution>(mFac.getStationaryProbs());
		_simulatedSequences = std::make_unique<sequenceContainer>();

        setRng(&_rng);

    }

	void setRng(RngType *rng) {
		_rng = rng;
	}

	void setSaveRates(bool saveRates) {
		_saveRates = saveRates;
	}


	void clearRatesVec() { 
		_siteRates.clear();
	}

	tree* gettree() {
		return _et;
	}

	std::unique_ptr<sequenceContainer> getSequenceContainer() {
		auto outputSequences = std::move(_simulatedSequences);
		_simulatedSequences = std::make_unique<sequenceContainer>();
		return outputSequences;
	}


	std::vector<double> getSiteRates() { 
		return _siteRates;
	}

	void generateSubstitutionsAlongTree(int seqLength) {
		std::vector<MDOUBLE> ratesVec(seqLength);
		MDOUBLE sumOfRatesAcrossSites = 0.0;
		_rateCategories.resize(seqLength);
		for (int h = 0; h < seqLength; h++)  {
			int selectedRandomCategory = _rateCategorySampler.drawSample(*_rng);
			_rateCategories[h] = selectedRandomCategory;
			ratesVec[h] = _sp->rates(selectedRandomCategory);
			sumOfRatesAcrossSites += ratesVec[h];
		}
		_siteRates.clear();
		if (_saveRates) _siteRates.insert(_siteRates.end(), ratesVec.begin(), ratesVec.end());

		sequence rootSequence = generateRootSeq(seqLength, ratesVec);

		if ((*_nodesToSave)[_et->getRoot()->id()]){ 
			saveSequence(rootSequence);
		}

		mutateSeqRecuresively(rootSequence, _et->getRoot());
	}

	void mutateSeqRecuresively(const sequence& currentSequence, tree::nodeP currentNode) {
		if (currentNode->isLeaf()) return;

		for (auto &node: currentNode->getSons()) {
			sequence childSeq(currentSequence);
			childSeq.setID(node->id());
			childSeq.setName(node->name());
			mutateSeqAlongBranch(childSeq, node->dis2father());
			if ((*_nodesToSave)[node->id()]) saveSequence(childSeq);
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

    std::shared_ptr<sequenceContainer> simulateSubstitutions(size_t sequenceLength) {
        generateSubstitutionsAlongTree(sequenceLength);

        return getSequenceContainer();
    }

    void setAlignedSequenceMap(const MSA& msa) {
        const auto& alignedSeq = msa.getAlignedSequence();
        setAlignedSequenceMap(alignedSeq);
    }


private:

	sequence generateRootSeq(int seqLength, std::vector<MDOUBLE>& ratesVec) {
		sequence rootSeq(_alph);

		rootSeq.resize(seqLength);

		size_t rootID = _et->getRoot()->id();
		for (int i = 0; i < seqLength; i++) {
			ALPHACHAR newChar = _frequencySampler->drawSample(*_rng) - 1;
			rootSeq[i] = newChar;
		}
		rootSeq.setName(_et->getRoot()->name());
		rootSeq.setID(_et->getRoot()->id());
		return rootSeq;

	}

	void mutateSeqAlongBranch(sequence& currentSequence, const MDOUBLE& distToFather) {
		mutateEntireSeq(currentSequence, distToFather);
	}

	void mutateEntireSeq(sequence& currentSequence, const MDOUBLE& branchLength) {
		const int nodeId = currentSequence.id();
		BranchTransitionProbabilities<AlphabetSize> cachedPijt(branchLength, *_sp);
		
		// Check if this is a leaf we're saving (low memory mode)
		if (_alignedSequenceMap != nullptr && (*_nodesToSave)[nodeId]) {
			const std::vector<int>& gapStructure = _alignedSequenceMap->at(nodeId);
			
			size_t site = 0;
			for (int blockSize : gapStructure) {
				if (blockSize < 0) {
					// Gap block - skip these sites
					site += (-blockSize);
					continue;
				}
				// Non-gap block - mutate these sites
				for (int i = 0; i < blockSize; ++i, ++site) {
					ALPHACHAR parentChar = currentSequence[site];
					auto &Pijt = cachedPijt.getDistribution(_rateCategories[site], parentChar);
					ALPHACHAR nextChar = Pijt.drawSample(*_rng) - 1;
					currentSequence[site] = nextChar;
				}
			}
		} else {
			// Normal mode - mutate all sites
			for (size_t site = 0; site < currentSequence.seqLen(); ++site) {
				ALPHACHAR parentChar = currentSequence[site];
				auto &Pijt = cachedPijt.getDistribution(_rateCategories[site], parentChar);
				ALPHACHAR nextChar = Pijt.drawSample(*_rng) - 1;
				currentSequence[site] = nextChar;
			}
		}
	}


	void saveSequence(const sequence &currentSequence) {
		if (_finalMsaPath.size() > 0) {
			saveSequenceToDisk(currentSequence);
			return;
		}
		sequence temp(currentSequence);
		_simulatedSequences->add(temp);
	}

	void saveSequenceToDisk(const sequence &currentSequence) {
		const int nodeId = currentSequence.id();
		
		
		_outputFile << ">" << currentSequence.name() << "\n";
		
		// Get gap structure for this sequence
		if (_alignedSequenceMap != nullptr) {
			const std::vector<int>& gapStructure = _alignedSequenceMap->at(nodeId);
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


	void setSaveStateLeaves(const tree::nodeP &node) {
		for(auto &node: node->getSons()) {
			if (node->isLeaf()) (*_nodesToSave)[node->id()] = true;
			setSaveStateLeaves(node);
		}
	}


	tree* _et;
	std::shared_ptr<const stochasticProcess> _sp;
	const alphabet* _alph;

	std::shared_ptr<std::vector<bool>> _nodesToSave;
	bool _saveRates;
	std::vector<std::unique_ptr<DiscreteDistribution>> _gillespieSampler;

	std::vector<size_t> _rateCategories;
	std::vector<double> _siteRates;
	std::unique_ptr<sequenceContainer> _simulatedSequences;
	std::unique_ptr<DiscreteDistribution> _frequencySampler;

	CategorySampler _rateCategorySampler;
	std::string _finalMsaPath;

	std::array<std::string, AlphabetSize> _charLookup;

	const std::unordered_map<size_t, std::vector<int>>* _alignedSequenceMap = nullptr;

	RngType *_rng;
	std::ofstream _outputFile;

};

#endif