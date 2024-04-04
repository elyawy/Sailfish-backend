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

rateMatrixSim::rateMatrixSim(modelFactory& mFac) {
		_et = mFac.getTree();
		_sp = mFac.getStochasticProcess();
		_alphaSize = _sp->alphabetSize();

		_cpijGam = computePijGam();
		_cpijGam.fillPij(*_et, *_sp);


		_subManager = std::make_shared<substitutionManager>(_et->getNodesNum());

		_alph = mFac.getAlphabet();
		_rootSequence = std::make_shared<sequence>(_alph);
		
		_avgSubtitutionsPerSite = 0.0;
		_subtitutionsRatePerSite = 0.0;

		_rateCategoriesSum = 0.0;
		std::vector<MDOUBLE> rateProbs;
		for (int j = 0 ; j < _sp->categories(); ++j) {
			_rateCategoriesSum += _sp->rates(j);
			rateProbs.push_back(_sp->ratesProb(j));
		}
		for (int j = 0 ; j < _sp->categories(); ++j) {
			_normalizedRateCategories.push_back(_sp->rates(j)/_rateCategoriesSum);
		}
		_rateSampler = std::make_unique<DiscreteDistribution>(rateProbs);

		std::vector<MDOUBLE> frequencies;
		for (int j = 0; j < _alphaSize; ++j) {
			frequencies.push_back(_sp->freq(j));
		}
		_frequencySampler = std::make_unique<DiscreteDistribution>(frequencies);

	};
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

void rateMatrixSim::setSeed(size_t seed) {
	_seed = seed;
	_mt_rand.seed(seed);
}


void rateMatrixSim::generate_substitution_log(int seqLength) {
	using nodeP = tree::nodeP;
	// generateRootLog(seqLength);
	// _rootSequence = sequence(_alph);

	std::vector<MDOUBLE> rateNormalizedVec(seqLength);

	_rateVec.resize(seqLength);
	for (int h = 0; h < seqLength; h++)  {
		int theRanCat = _rateSampler->drawSample() - 1;
		_rateVec[h] = _sp->rates(theRanCat);
		rateNormalizedVec[h] = _normalizedRateCategories[theRanCat];
	}

	_siteSampler = std::make_unique<DiscreteDistribution>(rateNormalizedVec);
	// _siteSampler->setSeed(_seed);

	_rootSequence->resize(seqLength);
	generateRootSeq(seqLength);

	// std::cout << "starting subs simulation\n";
	// std::cin.get();

	mutateSeqRecuresively(_et->getRoot(), seqLength);

	// _subManager->printSubManager();
}

void rateMatrixSim::mutateSeqRecuresively(tree::nodeP currentNode, int seqLength) {

	for (auto &node: currentNode->getSons()) {
		mutateSeqAlongBranch(node, seqLength);
		mutateSeqRecuresively(node, seqLength);
		if (node->isLeaf()) saveSequence(node->id(), node->name());

		// std::cout << "Node: " << currentNode->id() << "\n";
		if (!_subManager->isEmpty(currentNode->id())) {
			_subManager->undoSubs(currentNode->id(), *_rootSequence, _rateVec, _sp.get());
		}
	}
}

void rateMatrixSim::mutateSeqAlongBranch(tree::nodeP currentNode, int seqLength) {
	const MDOUBLE distToFather = currentNode->dis2father();
	if (distToFather > 0.2) {
		mutateEntireSeq(currentNode, seqLength);
	} else {
		mutateSeqGillespie(currentNode, seqLength, distToFather);
	}
}


void rateMatrixSim::mutateEntireSeq(tree::nodeP currentNode, int seqLength) {
	const int nodeId = currentNode->id();
	const int parentId = currentNode->father()->id();
	for (size_t site = 0; site < seqLength; ++site) {
		ALPHACHAR parentChar = (*_subManager).getCharacter(parentId, site, *_rootSequence);
		ALPHACHAR nextChar = _cpijGam.getRandomChar(_rateVec[site], nodeId, parentChar);
		if (nextChar != parentChar){
			_subManager->handleEvent(parentId, site, parentChar, _rateVec, _sp.get(), *_rootSequence);
			if (currentNode->isLeaf()) {
				_subManager->handleEvent(nodeId, site, nextChar, _rateVec, _sp.get(), *_rootSequence);
			}
		}
	}
}


void rateMatrixSim::mutateSeqGillespie(tree::nodeP currentNode, int seqLength, MDOUBLE distToParent) {
	// std::cout << "mutating sequence using Gillespie!\n";

	const int nodeId = currentNode->id();
	const int parentId = currentNode->father()->id();
	MDOUBLE branchLength = distToParent;

	double lambdaParam = _subManager->getReactantsSum();
	std::exponential_distribution<double> distribution(lambdaParam);
	double waitingTime = distribution(_mt_rand);
	if (waitingTime < 0) {
		std::cout << branchLength << " " << lambdaParam << " " << waitingTime << "\n";
		errorMsg::reportError("waiting time is negative :(");
	}
	while (waitingTime < branchLength) {
		if (waitingTime < 0) errorMsg::reportError("waiting time is negative :(");

		int mutatedSite = _siteSampler->drawSample() - 1;
		ALPHACHAR parentChar = (*_subManager).getCharacter(parentId, mutatedSite, *_rootSequence);
		ALPHACHAR nextChar = _cpijGam.getRandomChar(_rateVec[mutatedSite], nodeId, parentChar);
		_subManager->handleEvent(parentId, mutatedSite, parentChar, _rateVec, _sp.get(), *_rootSequence);
		if (currentNode->isLeaf()) {
			_subManager->handleEvent(nodeId, mutatedSite, nextChar, _rateVec, _sp.get(), *_rootSequence);
		}

		lambdaParam = _subManager->getReactantsSum();
		branchLength = branchLength - waitingTime;
		std::exponential_distribution<double> distribution(lambdaParam);
		waitingTime = distribution(_mt_rand);

	}
}


void rateMatrixSim::generateRootSeq(int seqLength) {	
	for (int i = 0; i < seqLength; i++) {
		ALPHACHAR newChar = _frequencySampler->drawSample() - 1;
		(*_rootSequence)[i] =  newChar;
		MDOUBLE qii = _sp->Qij(newChar, newChar);
		MDOUBLE rate = _rateVec[i];
		if(qii > 0) errorMsg::reportError("Qii is positive!");
		if(rate < 0) errorMsg::reportError("rate is negative!");
		_subManager->updateReactantsSum((_sp->Qij(newChar, newChar)),_rateVec[i]);
     }

	_rootSequence->setAlphabet(_alph);
	_rootSequence->setName(_et->getRoot()->name());
	_rootSequence->setID(_et->getRoot()->id());
}


void rateMatrixSim::saveSequence(int nodeId, std::string name) {
	std::unique_ptr<sequence> temp = std::make_unique<sequence>(*_rootSequence);
	temp->setName(name);
	temp->setID(nodeId);
	_simulatedSequences.push_back(std::move(temp));
}

// sequenceContainer rateMatrixSim::toSeqData() {
// 	sequenceContainer myseqData;
// 	for (int i=0; i < _simulatedSequences.size(); ++i) {
// 		myseqData.add(*_simulatedSequences[i]);
// 	}
// 	return myseqData;
// }



std::unique_ptr<sequenceContainer> rateMatrixSim::toSeqDataWithoutInternalNodes() {
	std::unique_ptr<sequenceContainer> myseqData = std::make_unique<sequenceContainer>();
	// sequenceContainer myseqData;
	for (int i=0; i < _simulatedSequences.size(); ++i) {
		tree::nodeP theCurNode = _et->findNodeByName(_simulatedSequences[i]->name());
		if (theCurNode == NULL)
			errorMsg::reportError("could not find the specified name: " + _simulatedSequences[i]->name());
		if (theCurNode->isInternal()) continue;

		myseqData->add(*std::move(_simulatedSequences[i]));
	}
	return myseqData;
}
