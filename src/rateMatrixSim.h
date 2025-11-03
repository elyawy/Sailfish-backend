// $Id: simulateTree.h 8507 2010-08-12 15:20:59Z rubi $

#ifndef ___RATE_MATRIX_SIM
#define ___RATE_MATRIX_SIM

#include "../libs/Phylolib/includes/definitions.h"
#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/stochasticProcess.h"
#include "../libs/Phylolib/includes/sequenceContainer.h"

#include "modelFactory.h"
// #include "substitutionManager.h"
#include "CategorySampler.h"
#include "CachedTransitionProbabilities.h"

template<typename RngType = std::mt19937_64>
class rateMatrixSim {
public:
	explicit rateMatrixSim(modelFactory& mFac, std::shared_ptr<std::vector<bool>> nodesToSave) : 
		_et(mFac.getTree()), _sp(mFac.getStochasticProcess()), _alph(mFac.getAlphabet()), 
		_invariantSitesProportion(mFac.getInvariantSitesProportion()),
		_siteRateCorrelation(mFac.getSiteRateCorrelation()),
		_cachedPijt(*mFac.getTree(), *mFac.getStochasticProcess()),
		_nodesToSave(nodesToSave), _saveRates(false),
		_rateCategorySampler(buildRateCategoryProbs(mFac), mFac.getSiteRateCorrelation()),
		_writeFolder("") {
		
		size_t alphaSize = _sp->alphabetSize();
		
		std::vector<MDOUBLE> frequencies;
		for (int j = 0; j < alphaSize; ++j) {
			frequencies.push_back(_sp->freq(j));
		}
		_frequencySampler = std::make_unique<DiscreteDistribution>(frequencies);
		_simulatedSequences = std::make_unique<sequenceContainer>();
	}

	virtual ~rateMatrixSim() {}

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

	void generate_substitution_log(int seqLength) {
		std::vector<MDOUBLE> ratesVec(seqLength);
		MDOUBLE sumOfRatesAcrossSites = 0.0;
		_rateCategories.resize(seqLength);
		for (int h = 0; h < seqLength; h++)  {
			int selectedRandomCategory = _rateCategorySampler.drawSample(*_rng);
			_rateCategories[h] = selectedRandomCategory;
			if (selectedRandomCategory >= _sp->categories()) {
				ratesVec[h] = 0.0;
				continue;
			}
			ratesVec[h] = _sp->rates(selectedRandomCategory);
			sumOfRatesAcrossSites += ratesVec[h];
		}

		if (_saveRates) _siteRates.insert(_siteRates.end(), ratesVec.begin(), ratesVec.end());

		// _currentSequence->resize(seqLength);
		sequence rootSequence = generateRootSeq(seqLength, ratesVec);

		if ((*_nodesToSave)[_et->getRoot()->id()]){ 
			saveSequence(rootSequence);
		}

		mutateSeqRecuresively(rootSequence, _et->getRoot());

		// _subManager.clear();
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
			// if (!(*_nodesToSave)[node->id()]) _subManager.freeSequence();
			// if (!_subManager.isEmpty(node->id())) {
			// 	_subManager.undoSubs(node->id(), *_currentSequence, _rateCategories, _sp.get());
			// }
		}
	}

	void setWriteFolder(const std::string &folder) {
		_writeFolder = folder;
	}

private:
    static std::vector<MDOUBLE> buildRateCategoryProbs(modelFactory& mFac) {
        std::vector<MDOUBLE> rateCategoriesProbs;
        auto sp = mFac.getStochasticProcess();
        MDOUBLE invariantProp = mFac.getInvariantSitesProportion();
        
        for (int j = 0; j < sp->categories(); ++j) {
            MDOUBLE currentRateProb = sp->ratesProb(j);
            currentRateProb = currentRateProb * (1.0 - invariantProp);
            rateCategoriesProbs.push_back(currentRateProb);
        }
        if (invariantProp > 0.0) rateCategoriesProbs.push_back(invariantProp);
        
        return rateCategoriesProbs;
    }

	sequence generateRootSeq(int seqLength, std::vector<MDOUBLE>& ratesVec) {
		sequence rootSeq(_alph);

		rootSeq.resize(seqLength);

		size_t rootID = _et->getRoot()->id();
		for (int i = 0; i < seqLength; i++) {
			ALPHACHAR newChar = _frequencySampler->drawSample(*_rng) - 1;
			rootSeq[i] = newChar;
		}
		// _subManager.setRootSequence(seqLength, ratesVec, _sp.get(), *_currentSequence);
		
		// rootSeq.setAlphabet(_alph);
		rootSeq.setName(_et->getRoot()->name());
		rootSeq.setID(_et->getRoot()->id());
		return rootSeq;

	}

	void mutateSeqAlongBranch(sequence& currentSequence, const MDOUBLE& distToFather) {
		// const MDOUBLE distToFather = currentNode->dis2father();
		mutateEntireSeq(currentSequence);
	}

	void mutateEntireSeq(sequence& currentSequence) {
		const int nodeId = currentSequence.id();
		// const int parentId = currentNode->father()->id();

		for (size_t site = 0; site < currentSequence.seqLen(); ++site) {
			ALPHACHAR parentChar = currentSequence[site];
			if (_rateCategories[site] == _sp->categories()) continue;
			auto &Pijt = _cachedPijt.getDistribution(nodeId, _rateCategories[site], parentChar);
			ALPHACHAR nextChar = Pijt.drawSample(*_rng) - 1;
			currentSequence[site] = nextChar;
			// if (nextChar != parentChar){
			// 	_subManager.handleEvent(nodeId, site, nextChar, _rateCategories, _sp.get(), *_currentSequence);
			// }
		}
	}

	// void mutateSeqGillespie(tree::nodeP currentNode, int seqLength, MDOUBLE distToParent) {
	// 	const int nodeId = currentNode->id();
	// 	const int parentId = currentNode->father()->id();
	// 	MDOUBLE branchLength = distToParent;

	// 	double lambdaParam = _subManager.getReactantsSum();
	// 	std::exponential_distribution<double> distribution(-lambdaParam);
	// 	double waitingTime = distribution(*_rng);
	// 	if (waitingTime < 0) {
	// 		std::cout << branchLength << " " << lambdaParam << " " << waitingTime << "\n";
	// 		errorMsg::reportError("waiting time is negative :(");
	// 	}
	// 	while (waitingTime < branchLength) {
	// 		if (waitingTime < 0) {
	// 			std::cout << branchLength << " " << lambdaParam << " " << waitingTime << "\n";
	// 			errorMsg::reportError("waiting time is negative :(");
	// 		}

	// 		int mutatedSite = _subManager.sampleSite(*_rng);
	// 		ALPHACHAR parentChar = (*_currentSequence)[mutatedSite];
	// 		ALPHACHAR nextChar = _gillespieSampler[parentChar]->drawSample(*_rng) - 1;
	// 		_subManager.handleEvent(nodeId, mutatedSite, nextChar, _rateCategories, _sp.get(), *_currentSequence);

	// 		lambdaParam = _subManager.getReactantsSum();
	// 		branchLength = branchLength - waitingTime;
	// 		std::exponential_distribution<double> distribution(-lambdaParam);
	// 		waitingTime = distribution(*_rng);
	// 	}
	// }

	void saveSequence(const sequence &currentSequence) {
		if (_writeFolder.size() > 0) {
			saveSequenceToDisk(currentSequence);
			return;
		}
		sequence temp(currentSequence);
		_simulatedSequences->add(temp);
	}

	void saveSequenceToDisk(const sequence &currentSequence) {
		std::string fileName = _writeFolder + "/" + std::to_string(currentSequence.id()) + ".fasta";
		std::ofstream outFile(fileName);
		if (!outFile.is_open()) {
			errorMsg::reportError("Could not open file " + fileName + " for writing simulated sequences.");
		}
		outFile << ">" << currentSequence.name() << "\n";
		outFile << currentSequence.toString() << "\n";
		outFile.close();
	}

	void initGillespieSampler() {
		_gillespieSampler.resize(_alph->size());
		for (size_t i = 0; i < _alph->size(); ++i) {
			std::vector<double> qRates(_alph->size(), 0.0);
			double sum = -_sp->Qij(i,i);
			double normalizer = 1.0 / sum;
			for (size_t j = 0; j < _alph->size(); ++j) {
				if (i==j) continue;
				qRates[j] = _sp->Qij(i,j) * normalizer;
			}
			_gillespieSampler[i] = std::make_unique<DiscreteDistribution>(qRates);
		}
	}

	void setSaveStateLeaves(const tree::nodeP &node) {
		for(auto &node: node->getSons()) {
			if (node->isLeaf()) (*_nodesToSave)[node->id()] = true;
			setSaveStateLeaves(node);
		}
	}

	// bool testSumOfRates() {
	// 	MDOUBLE sumOfRates = 0.0;
	// 	for (size_t i = 0; i < _currentSequence->seqLen(); i++) {
	// 		ALPHACHAR currentChar = (*_currentSequence)[i];
	// 		MDOUBLE currentQii = _sp->Qij(currentChar, currentChar);
	// 		MDOUBLE currentRate = _sp->rates(_rateCategories[i]);
	// 		sumOfRates += (currentQii*currentRate);
	// 	}
	// 	MDOUBLE preCalculatedSum = _subManager.getReactantsSum();
	// 	if (abs(preCalculatedSum - sumOfRates) > 1e-6) {
	// 		std::cout << "preCalculatedSum=" << preCalculatedSum << " "
	// 				  << "sumOfRates=" << sumOfRates;
	// 		errorMsg::reportError("Error in sum of rates calculation!");
	// 	}
	// 	std::cout << "preCalculatedSum is correct\n" << "preCalculatedSum=" << preCalculatedSum << " "
	// 				  << "sumOfRates=" << sumOfRates << "\n";

	// 	return true;
	// }

	tree* _et;
	std::shared_ptr<const stochasticProcess> _sp;
	const alphabet* _alph;

	MDOUBLE _invariantSitesProportion;
	MDOUBLE _siteRateCorrelation;

	CachedTransitionProbabilities _cachedPijt;
	// computePijGam _cpijGam;
	// sequence* _currentSequence;
	// substitutionManager _subManager;
	std::shared_ptr<std::vector<bool>> _nodesToSave;
	bool _saveRates;
	std::vector<std::unique_ptr<DiscreteDistribution>> _gillespieSampler;

	std::vector<size_t> _rateCategories;
	std::vector<double> _siteRates;
	std::unique_ptr<sequenceContainer> _simulatedSequences;
	std::unique_ptr<DiscreteDistribution> _frequencySampler;

	CategorySampler _rateCategorySampler;
	std::string _writeFolder;

	RngType *_rng;
};

#endif