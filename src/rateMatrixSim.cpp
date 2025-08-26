// $Id: simulateTree.cpp 8508 2010-08-12 15:21:04Z rubi $
#include <stack>
#include <unordered_map>
#include <ostream>
#include <sstream>


#include "../libs/Phylolib/includes/definitions.h"
#include "../libs/Phylolib/includes/treeUtil.h"
#include "../libs/Phylolib/includes/talRandom.h"
#include "../libs/Phylolib/includes/gammaDistribution.h"
#include "../libs/Phylolib/includes/codon.h"

#include "rateMatrixSim.h"
// simulateTree::simulateTree(tree*  _inEt,
// 						   const stochasticProcess* sp,
// 						   const alphabet* alph) :
// 	_et(_inEt), _sp(sp),_alph(alph),_avgSubtitutionsPerSite(0.0) {
// 	};

rateMatrixSim::rateMatrixSim(modelFactory& mFac, std::shared_ptr<std::vector<bool>> nodesToSave) : 
	_et(mFac.getTree()), _sp(mFac.getStochasticProcess()), _alph(mFac.getAlphabet()), 
	_invariantSitesProportion(mFac.getInvariantSitesProportion()),
	_cpijGam(), _rootSequence(mFac.getAlphabet()), _subManager(mFac.getTree()->getNodesNum()),
	_nodesToSave(nodesToSave), _saveRates(false), _biased_coin(0,1) {
		// _et = mFac.getTree();
		// _sp = mFac.getStochasticProcess();
		// _alph = mFac.getAlphabet();

		size_t alphaSize = _sp->alphabetSize();

		_cpijGam.fillPij(*_et, *_sp);
		initGillespieSampler();
		

		std::vector<MDOUBLE> rateProbs;
		for (int j = 0 ; j < _sp->categories(); ++j) {
			MDOUBLE currentRateProb = _sp->ratesProb(j);
			currentRateProb = currentRateProb * (1.0  - _invariantSitesProportion);
			rateProbs.push_back(currentRateProb);
		}
		if (_invariantSitesProportion > 0.0) rateProbs.push_back(_invariantSitesProportion);

		_rateSampler = std::make_unique<DiscreteDistribution>(rateProbs);

		std::vector<MDOUBLE> frequencies;
		for (int j = 0; j < alphaSize; ++j) {
			frequencies.push_back(_sp->freq(j));
		}
		_frequencySampler = std::make_unique<DiscreteDistribution>(frequencies);

		_simulatedSequences = std::make_unique<sequenceContainer>();

};

void rateMatrixSim::setSaveRates(bool saveRates) {
	_saveRates = saveRates;
}

void rateMatrixSim::initGillespieSampler() {
	_gillespieSampler.resize(_alph->size());
	for (size_t i = 0; i < _alph->size(); ++i) {
		std::vector<double> qRates(_alph->size(), 0.0);
		double sum = -_sp->Qij(i,i);
		double normalizer = 1.0 / sum;
		for (size_t j = 0; j < _alph->size(); ++j) {
			if (i==j) continue;
			qRates[j] = _sp->Qij(i,j) * normalizer;
			// std::cout << i << j << "->" << qRates[j] << ",";
		}
		// std::cout << "\n" << i << " " << sum << "\n";
		_gillespieSampler[i] = std::make_unique<DiscreteDistribution>(qRates);
	}
}

// simulateTree::simulateTree(const tree&  _inEt,
// 						   const stochasticProcess& sp,
// 						   const alphabet* alph) : _sp(sp) {
// 		_et = _inEt;
// 		// _sp = sp;
// 		_alph = alph;
// 		_avgSubtitutionsPerSite = 0.0;
// 	};

rateMatrixSim::~rateMatrixSim() {
}

// void rateMatrixSim::setSeed(size_t seed) {
// 	_seed = seed;
// 	_mt_rand->seed(seed);
// }

void rateMatrixSim::setRng(mt19937_64 *rng) {
	_mt_rand = rng;
}


// const mt19937_64& rateMatrixSim::getRng(){
// 	return *_mt_rand;
// }


void rateMatrixSim::generate_substitution_log(int seqLength) {
	std::vector<MDOUBLE> ratesVec(seqLength);

	MDOUBLE sumOfRatesAcrossSites = 0.0;
	_rateCategories.resize(seqLength);
	for (int h = 0; h < seqLength; h++)  {
		int selectedRandomCategory = _rateSampler->drawSample() - 1;
		_rateCategories[h] = selectedRandomCategory;
		if (selectedRandomCategory >= _sp->categories()) {
			ratesVec[h] = 0.0;
			continue;
		}
		ratesVec[h] = _sp->rates(selectedRandomCategory);
		sumOfRatesAcrossSites += ratesVec[h];
	}
	if (_saveRates) _siteRates.insert(_siteRates.end(), ratesVec.begin(), ratesVec.end());
	// MDOUBLE sumOfRatesNoramlizingFactor = 1.0 / sumOfRatesAcrossSites;

	// _siteSampler = std::make_unique<DiscreteDistribution>(ratesVec, sumOfRatesNoramlizingFactor);
	_rootSequence.resize(seqLength);
	generateRootSeq(seqLength, ratesVec);
	if ((*_nodesToSave)[_et->getRoot()->id()]) saveSequence(_et->getRoot()->id(), _et->getRoot()->name());

	mutateSeqRecuresively(_et->getRoot(), seqLength);
	_subManager.clear();
}

void rateMatrixSim::mutateSeqRecuresively(tree::nodeP currentNode, int seqLength) {
	if (currentNode->isLeaf()) return;

	for (auto &node: currentNode->getSons()) {
		mutateSeqAlongBranch(node, seqLength);
		if ((*_nodesToSave)[node->id()]) saveSequence(node->id(), node->name());
		mutateSeqRecuresively(node, seqLength);
		if (!_subManager.isEmpty(node->id())) {
			_subManager.undoSubs(node->id(), _rootSequence, _rateCategories, _sp.get());
		}
	}
}

void rateMatrixSim::mutateSeqAlongBranch(tree::nodeP currentNode, int seqLength) {
	const MDOUBLE distToFather = currentNode->dis2father();
	mutateEntireSeq(currentNode, seqLength);
	
	// if (distToFather > 0.5) {
	// 	mutateEntireSeq(currentNode, seqLength);
	// } else {
	// 	mutateSeqGillespie(currentNode, seqLength, distToFather);
	// }
	// testSumOfRates();
}


void rateMatrixSim::mutateEntireSeq(tree::nodeP currentNode, int seqLength) {
	const int nodeId = currentNode->id();
	const int parentId = currentNode->father()->id();

	for (size_t site = 0; site < seqLength; ++site) {
		ALPHACHAR parentChar = _rootSequence[site];//_subManager.getCharacter(parentId, site, _rootSequence);
		if (_rateCategories[site] == _sp->categories()) continue;
		ALPHACHAR nextChar = _cpijGam.getRandomChar(_rateCategories[site], nodeId, parentChar);
		if (nextChar != parentChar){
			_subManager.handleEvent(nodeId, site, nextChar, _rateCategories, _sp.get(), _rootSequence);
		}
	}
}


void rateMatrixSim::mutateSeqGillespie(tree::nodeP currentNode, int seqLength, MDOUBLE distToParent) {
	// std::cout << "mutating sequence using Gillespie!\n";

	const int nodeId = currentNode->id();
	const int parentId = currentNode->father()->id();
	MDOUBLE branchLength = distToParent;

	double lambdaParam = _subManager.getReactantsSum();
	std::exponential_distribution<double> distribution(-lambdaParam);
	double waitingTime = distribution(*_mt_rand);
	if (waitingTime < 0) {
		std::cout << branchLength << " " << lambdaParam << " " << waitingTime << "\n";
		errorMsg::reportError("waiting time is negative :(");
	}
	while (waitingTime < branchLength) {
		if (waitingTime < 0) {
			std::cout << branchLength << " " << lambdaParam << " " << waitingTime << "\n";
			errorMsg::reportError("waiting time is negative :(");
		}

		int mutatedSite = _subManager.sampleSite(*_mt_rand);
		ALPHACHAR parentChar = _rootSequence[mutatedSite];
		ALPHACHAR nextChar = _gillespieSampler[parentChar]->drawSample() - 1;
		// std::cout << (int)parentChar << "->" << (int)nextChar << "\n";
		_subManager.handleEvent(nodeId, mutatedSite, nextChar, _rateCategories, _sp.get(), _rootSequence);

		lambdaParam = _subManager.getReactantsSum();
		branchLength = branchLength - waitingTime;
		std::exponential_distribution<double> distribution(-lambdaParam);
		waitingTime = distribution(*_mt_rand);

	}
}




void rateMatrixSim::generateRootSeq(int seqLength, std::vector<MDOUBLE>& ratesVec) {
	size_t rootID = _et->getRoot()->id();
	for (int i = 0; i < seqLength; i++) {
		ALPHACHAR newChar = _frequencySampler->drawSample() - 1;
		// ratesVec[i] = ratesVec[i]*(-_sp->Qij(newChar, newChar));
		_rootSequence[i] =  newChar;
     }
	// std::cout << ">Root-sequence\n" << _rootSequence  <<  "\n";
	// std::cout << ">Rates\n" << ratesVec;
	_subManager.handleRootSequence(seqLength, ratesVec, _sp.get(), _rootSequence);
	
	_rootSequence.setAlphabet(_alph);
	_rootSequence.setName(_et->getRoot()->name());
	_rootSequence.setID(_et->getRoot()->id());
}


void rateMatrixSim::saveSequence(const int &nodeId,const std::string &name) {
	sequence temp(_rootSequence);
	temp.setName(name);
	temp.setID(nodeId);
	// std::cout << temp << "\n";
	_simulatedSequences->add(temp);
}

// sequenceContainer rateMatrixSim::toSeqData() {
// 	sequenceContainer myseqData;
// 	for (int i=0; i < _simulatedSequences.size(); ++i) {
// 		myseqData.add(*_simulatedSequences[i]);
// 	}
// 	return myseqData;
// }



std::unique_ptr<sequenceContainer> rateMatrixSim::getSequenceContainer() {
	// std::unique_ptr<sequenceContainer> myseqData = std::make_unique<sequenceContainer>();
	// // sequenceContainer myseqData;
	// for (int i=0; i < _simulatedSequences.size(); ++i) {
	// 	tree::nodeP theCurNode = _et->findNodeById(_simulatedSequences[i]->id());
	// 	if (theCurNode == NULL)
	// 		errorMsg::reportError("could not find the specified name: " + _simulatedSequences[i]->name());
	// 	if (theCurNode->isInternal()) continue;
	auto outputSequences = std::move(_simulatedSequences);
	_simulatedSequences = std::make_unique<sequenceContainer>();
	// 	myseqData->add(*std::move(_simulatedSequences[i]));
	// }

	return outputSequences;
}


bool rateMatrixSim::testSumOfRates() {
	MDOUBLE sumOfRates = 0.0;
	for (size_t i = 0; i < _rootSequence.seqLen(); i++) {
		ALPHACHAR currentChar = _rootSequence[i];
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
