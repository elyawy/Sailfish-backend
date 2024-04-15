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

rateMatrixSim::rateMatrixSim(modelFactory& mFac) : 
	_et(mFac.getTree()), _sp(mFac.getStochasticProcess()), _alph(mFac.getAlphabet()),
	_cpijGam(), _rootSequence(mFac.getAlphabet()), _subManager(mFac.getTree()->getNodesNum()),
	_nodesToSave(_et->getNodesNum(), false) {
		// _et = mFac.getTree();
		// _sp = mFac.getStochasticProcess();
		// _alph = mFac.getAlphabet();
		setSaveStateLeaves(_et->getRoot());

		size_t alphaSize = _sp->alphabetSize();

		// _cpijGam = computePijGam();
		_cpijGam.fillPij(*_et, *_sp);
		// _rootSequence = sequence(_alph);
		

		std::vector<MDOUBLE> rateProbs;
		for (int j = 0 ; j < _sp->categories(); ++j) {
			rateProbs.push_back(_sp->ratesProb(j));
		}
		_rateSampler = std::make_unique<DiscreteDistribution>(rateProbs);

		std::vector<MDOUBLE> frequencies;
		for (int j = 0; j < alphaSize; ++j) {
			frequencies.push_back(_sp->freq(j));
		}
		_frequencySampler = std::make_unique<DiscreteDistribution>(frequencies);

		_simulatedSequences = std::make_unique<sequenceContainer>();

};

void rateMatrixSim::setSaveStateLeaves(const tree::nodeP &node) {
	for(auto &node: node->getSons()) {
		if (node->isLeaf()) _nodesToSave[node->id()] = true;
		setSaveStateLeaves(node);
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
	using nodeP = tree::nodeP;
	// generateRootLog(seqLength);
	// _rootSequence = sequence(_alph);

	std::vector<MDOUBLE> ratesVec(seqLength);

	// _rateVec.resize(seqLength);
	MDOUBLE sumOfRatesAcrossSites = 0.0;
	_rateCategories.resize(seqLength);
	for (int h = 0; h < seqLength; h++)  {
		int selectedRandomCategory = _rateSampler->drawSample() - 1;
		_rateCategories[h] = selectedRandomCategory;
		ratesVec[h] = _sp->rates(selectedRandomCategory);
		sumOfRatesAcrossSites += ratesVec[h];
	}
	MDOUBLE sumOfRatesNoramlizingFactor = 1.0 / sumOfRatesAcrossSites;

	_siteSampler = std::make_unique<DiscreteDistribution>(ratesVec, sumOfRatesNoramlizingFactor);
	// _siteSampler->setSeed(_seed);

	_rootSequence.resize(seqLength);
	generateRootSeq(seqLength);
	// std::cout << "starting subs simulation\n";
	// std::cin.get();
	if (_nodesToSave[_et->getRoot()->id()]) saveSequence(_et->getRoot()->id(), _et->getRoot()->name());

	mutateSeqRecuresively(_et->getRoot(), seqLength);
	_subManager.clear();
	// _subManager.printSubManager();
}

void rateMatrixSim::mutateSeqRecuresively(tree::nodeP currentNode, int seqLength) {

	for (auto &node: currentNode->getSons()) {
		mutateSeqAlongBranch(node, seqLength);
		mutateSeqRecuresively(node, seqLength);
		if (_nodesToSave[node->id()]) saveSequence(node->id(), node->name());

		// std::cout << "Node: " << currentNode->id() << "\n";
		if (!_subManager.isEmpty(currentNode->id())) {
			_subManager.undoSubs(currentNode->id(), _rootSequence, _rateCategories, _sp.get());
		}
	}
}

void rateMatrixSim::mutateSeqAlongBranch(tree::nodeP currentNode, int seqLength) {
	const MDOUBLE distToFather = currentNode->dis2father();
	if (distToFather >10.0) {
		mutateEntireSeq(currentNode, seqLength);
	} else {
		mutateSeqGillespie(currentNode, seqLength, distToFather);
	}
}


void rateMatrixSim::mutateEntireSeq(tree::nodeP currentNode, int seqLength) {
	const int nodeId = currentNode->id();
	const int parentId = currentNode->father()->id();

	for (size_t site = 0; site < seqLength; ++site) {
		ALPHACHAR parentChar = _rootSequence[site];//_subManager.getCharacter(parentId, site, _rootSequence);
		ALPHACHAR nextChar = _cpijGam.getRandomChar(_rateCategories[site], nodeId, parentChar);
		if (nextChar != parentChar){
			_subManager.handleEvent(parentId, site, parentChar, _rateCategories, _sp.get(), _rootSequence);
			if (currentNode->isLeaf()) {
				_subManager.handleEvent(nodeId, site, nextChar, _rateCategories, _sp.get(), _rootSequence);
			}
		}
	}
}


void rateMatrixSim::mutateSeqGillespie(tree::nodeP currentNode, int seqLength, MDOUBLE distToParent) {
	// std::cout << "mutating sequence using Gillespie!\n";

	const int nodeId = currentNode->id();
	const int parentId = currentNode->father()->id();
	MDOUBLE branchLength = distToParent;

	double lambdaParam = _subManager.getReactantsSum();
	std::exponential_distribution<double> distribution(lambdaParam);
	double waitingTime = distribution(*_mt_rand);
	if (waitingTime < 0) {
		std::cout << branchLength << " " << lambdaParam << " " << waitingTime << "\n";
		errorMsg::reportError("waiting time is negative :(");
	}
	while (waitingTime < branchLength) {
		if (waitingTime < 0) errorMsg::reportError("waiting time is negative :(");

		int mutatedSite = _siteSampler->drawSample() - 1;
		ALPHACHAR parentChar = _rootSequence[mutatedSite];//_subManager.getCharacter(parentId, mutatedSite, _rootSequence);
		ALPHACHAR nextChar = _cpijGam.getRandomChar(_rateCategories[mutatedSite], nodeId, parentChar);
		_subManager.handleEvent(parentId, mutatedSite, parentChar, _rateCategories, _sp.get(), _rootSequence);
		if (currentNode->isLeaf()) {
			_subManager.handleEvent(nodeId, mutatedSite, nextChar, _rateCategories, _sp.get(), _rootSequence);
		}

		lambdaParam = _subManager.getReactantsSum();
		branchLength = branchLength - waitingTime;
		std::exponential_distribution<double> distribution(lambdaParam);
		waitingTime = distribution(*_mt_rand);

	}
}


void rateMatrixSim::generateRootSeq(int seqLength) {	
	for (int i = 0; i < seqLength; i++) {
		ALPHACHAR newChar = _frequencySampler->drawSample() - 1;
		_rootSequence[i] =  newChar;
		MDOUBLE qii = _sp->Qij(newChar, newChar);
		size_t rateCategory = _rateCategories[i];
		if(qii > 0) errorMsg::reportError("Qii is positive!");
		if(rateCategory < 0) errorMsg::reportError("rate category is negative!");
		_subManager.updateReactantsSum((_sp->Qij(newChar, newChar)),_sp->rates(rateCategory));
     }

	_rootSequence.setAlphabet(_alph);
	_rootSequence.setName(_et->getRoot()->name());
	_rootSequence.setID(_et->getRoot()->id());
}


void rateMatrixSim::saveSequence(const int &nodeId,const std::string &name) {
	std::unique_ptr<sequence> temp = std::make_unique<sequence>(_rootSequence);
	temp->setName(name);
	temp->setID(nodeId);
	_simulatedSequences->add(*std::move(temp));
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

	return std::move(outputSequences);
}

void rateMatrixSim::setNodesToSaves(std::vector<size_t> nodeIDs) {
	std::fill(_nodesToSave.begin(), _nodesToSave.end(), false);
	for(auto &nodeID: nodeIDs) {
		_nodesToSave[nodeID] = true;
	}
}
