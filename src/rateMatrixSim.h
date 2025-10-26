// $Id: simulateTree.h 8507 2010-08-12 15:20:59Z rubi $

#ifndef ___RATE_MATRIX_SIM
#define ___RATE_MATRIX_SIM

#include "../libs/Phylolib/includes/definitions.h"
#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/stochasticProcess.h"
#include "../libs/Phylolib/includes/sequenceContainer.h"
#include "../libs/Phylolib/includes/computePijComponent.h"

#include "modelFactory.h"
#include "substitutionManager.h"
#include "CategorySampler.h"

template<typename RngType = std::mt19937_64>
class rateMatrixSim {
public:
	explicit rateMatrixSim(modelFactory& mFac, std::shared_ptr<std::vector<bool>> nodesToSave) : 
		_et(mFac.getTree()), _sp(mFac.getStochasticProcess()), _alph(mFac.getAlphabet()), 
		_invariantSitesProportion(mFac.getInvariantSitesProportion()),
		_siteRateCorrelation(mFac.getSiteRateCorrelation()),
		_cpijGam(), _currentSequence(std::make_unique<sequence>(mFac.getAlphabet())), 
		_subManager(mFac.getTree()->getNodesNum()),
		_nodesToSave(nodesToSave), _saveRates(false),
		_rateCategorySampler(buildRateCategoryProbs(mFac), mFac.getSiteRateCorrelation()) {
		
		size_t alphaSize = _sp->alphabetSize();

		_cpijGam.fillPij(*_et, *_sp);
		initGillespieSampler();
				
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

		_currentSequence->resize(seqLength);
		generateRootSeq(seqLength, ratesVec);

		if ((*_nodesToSave)[_et->getRoot()->id()]) saveSequence(_et->getRoot()->id(), _et->getRoot()->name());

		mutateSeqRecuresively(_et->getRoot(), seqLength);

		_subManager.clear();
	}

	void mutateSeqRecuresively(tree::nodeP currentNode, int seqLength) {
		if (currentNode->isLeaf()) return;

		for (auto &node: currentNode->getSons()) {
			mutateSeqAlongBranch(node, seqLength);
			if ((*_nodesToSave)[node->id()]) saveSequence(node->id(), node->name());
			mutateSeqRecuresively(node, seqLength);
			if (!_subManager.isEmpty(node->id())) {
				_subManager.undoSubs(node->id(), *_currentSequence, _rateCategories, _sp.get());
			}
		}
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

	void generateRootSeq(int seqLength, std::vector<MDOUBLE>& ratesVec) {
		size_t rootID = _et->getRoot()->id();
		for (int i = 0; i < seqLength; i++) {
			ALPHACHAR newChar = _frequencySampler->drawSample(*_rng) - 1;
			(*_currentSequence)[i] = newChar;
		}
		_subManager.handleRootSequence(seqLength, ratesVec, _sp.get(), *_currentSequence);
		
		_currentSequence->setAlphabet(_alph);
		_currentSequence->setName(_et->getRoot()->name());
		_currentSequence->setID(_et->getRoot()->id());
	}

	void mutateSeqAlongBranch(tree::nodeP currentNode, int seqLength) {
		const MDOUBLE distToFather = currentNode->dis2father();
		mutateEntireSeq(currentNode, seqLength);
	}

	void mutateEntireSeq(tree::nodeP currentNode, int seqLength) {
		const int nodeId = currentNode->id();
		const int parentId = currentNode->father()->id();

		for (size_t site = 0; site < seqLength; ++site) {
			ALPHACHAR parentChar = (*_currentSequence)[site];
			if (_rateCategories[site] == _sp->categories()) continue;
			ALPHACHAR nextChar = _cpijGam.getRandomChar(_rateCategories[site], nodeId, parentChar, *_rng);
			if (nextChar != parentChar){
				_subManager.handleEvent(nodeId, site, nextChar, _rateCategories, _sp.get(), *_currentSequence);
			}
		}
	}

	void mutateSeqGillespie(tree::nodeP currentNode, int seqLength, MDOUBLE distToParent) {
		const int nodeId = currentNode->id();
		const int parentId = currentNode->father()->id();
		MDOUBLE branchLength = distToParent;

		double lambdaParam = _subManager.getReactantsSum();
		std::exponential_distribution<double> distribution(-lambdaParam);
		double waitingTime = distribution(*_rng);
		if (waitingTime < 0) {
			std::cout << branchLength << " " << lambdaParam << " " << waitingTime << "\n";
			errorMsg::reportError("waiting time is negative :(");
		}
		while (waitingTime < branchLength) {
			if (waitingTime < 0) {
				std::cout << branchLength << " " << lambdaParam << " " << waitingTime << "\n";
				errorMsg::reportError("waiting time is negative :(");
			}

			int mutatedSite = _subManager.sampleSite(*_rng);
			ALPHACHAR parentChar = (*_currentSequence)[mutatedSite];
			ALPHACHAR nextChar = _gillespieSampler[parentChar]->drawSample(*_rng) - 1;
			_subManager.handleEvent(nodeId, mutatedSite, nextChar, _rateCategories, _sp.get(), *_currentSequence);

			lambdaParam = _subManager.getReactantsSum();
			branchLength = branchLength - waitingTime;
			std::exponential_distribution<double> distribution(-lambdaParam);
			waitingTime = distribution(*_rng);
		}
	}

	void saveSequence(const int &nodeId, const std::string &name) {
		sequence temp(*_currentSequence);
		temp.setName(name);
		temp.setID(nodeId);
		_simulatedSequences->add(temp);
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

	bool testSumOfRates() {
		MDOUBLE sumOfRates = 0.0;
		for (size_t i = 0; i < _currentSequence->seqLen(); i++) {
			ALPHACHAR currentChar = (*_currentSequence)[i];
			MDOUBLE currentQii = _sp->Qij(currentChar, currentChar);
			MDOUBLE currentRate = _sp->rates(_rateCategories[i]);
			sumOfRates += (currentQii*currentRate);
		}
		MDOUBLE preCalculatedSum = _subManager.getReactantsSum();
		if (abs(preCalculatedSum - sumOfRates) > 1e-6) {
			std::cout << "preCalculatedSum=" << preCalculatedSum << " "
					  << "sumOfRates=" << sumOfRates;
			errorMsg::reportError("Error in sum of rates calculation!");
		}
		std::cout << "preCalculatedSum is correct\n" << "preCalculatedSum=" << preCalculatedSum << " "
					  << "sumOfRates=" << sumOfRates << "\n";

		return true;
	}

	tree* _et;
	std::shared_ptr<const stochasticProcess> _sp;
	const alphabet* _alph;

	MDOUBLE _invariantSitesProportion;
	MDOUBLE _siteRateCorrelation;

	computePijGam _cpijGam;
	std::unique_ptr<sequence> _currentSequence;
	substitutionManager _subManager;
	std::shared_ptr<std::vector<bool>> _nodesToSave;
	bool _saveRates;
	std::vector<std::unique_ptr<DiscreteDistribution>> _gillespieSampler;

	std::vector<size_t> _rateCategories;
	std::vector<double> _siteRates;
	std::unique_ptr<sequenceContainer> _simulatedSequences;
	std::unique_ptr<DiscreteDistribution> _frequencySampler;

	CategorySampler _rateCategorySampler;

	RngType *_rng;
};

#endif