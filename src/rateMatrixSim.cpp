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
		for (int j=0;j<_sp->categories() ;++j) {
			_rateCategoriesSum += _sp->rates(j);
		}
		for (int j=0;j<_sp->categories() ;++j) {
			_normalizedRateCategories.push_back(_sp->rates(j)/_rateCategoriesSum);
		}
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
	_cpijGam.setSeed(seed);
	talRandom::setSeed(seed);
}

// void rateMatrixSim::generate_seq(int seqLength) {

// 	sequence justAseq(_alph);
// 	_simulatedSequences.resize(_et->getNodesNum(),justAseq);

// 	for (int i=0; i < _simulatedSequences.size(); ++i) {
// 		_simulatedSequences[i].resize(seqLength);
// 	}

// 	generateRootSeq(seqLength);

// 	vector<MDOUBLE> rateVec(seqLength);
// 	for (int h = 0; h < seqLength; h++)  {
// 		int theRanCat = getRandCategory(h);
// 		rateVec[h] = _sp->rates(theRanCat);
// 	}
	
// 	_avgSubtitutionsPerSite = 0.0;
// 	for (int p=0 ; p < _et->getRoot()->getNumberOfSons() ; ++p) {
// 	  recursiveGenerateSpecificSeq(rateVec, seqLength, _et->getRoot()->getSon(p));
// 	}
// 	_avgSubtitutionsPerSite /= 1.0*seqLength;
// }

void rateMatrixSim::generate_substitution_log(int seqLength) {
	using nodeP = tree::nodeP;
	// generateRootLog(seqLength);
	// _rootSequence = sequence(_alph);

	std::vector<MDOUBLE> rateNormalizedVec(seqLength);

	_rateVec.resize(seqLength);
	for (int h = 0; h < seqLength; h++)  {
		int theRanCat = getRandCategory(h);
		_rateVec[h] = _sp->rates(theRanCat);
		rateNormalizedVec[h] = _normalizedRateCategories[theRanCat];
	}

	_siteSampler = std::make_unique<DiscreteDistribution>(rateNormalizedVec);

	_rootSequence->resize(seqLength);
	generateRootSeq(seqLength);

	std::cout << "starting subs simulation\n";
	std::cin.get();

	mutateSeqRecuresively(_et->getRoot(), seqLength);

	// _subManager->printSubManager();
}

void rateMatrixSim::mutateSeqRecuresively(tree::nodeP currentNode, int seqLength) {

	for (auto &node: currentNode->getSons()) {
		// std::cout << "mutating from " << currentNode->id() << " to " << node->id() << "\n";


		mutateSeqAlongBranch(node, seqLength);
		mutateSeqRecuresively(node, seqLength);
		// std::cout << "finished from " << currentNode->id() << " to " << node->id() << "\n";

		// if (node->isLeaf()) {
			// std::stringstream filename;
			// filename << "/home/elyawy/temp/" << node->id() << ".fasta";
			// std::cout << "writing to " << filename.str() << "\n";

			// std::ofstream o(filename.str());
			// o <<  (*_rootSequence) << std::endl;
			// o.close();
			// _subManager->dumpSubstitutionLog(node->id());
		// }
		std::cout << currentNode->id() << "\n";
		// _subManager->dumpSubstitutionLog(currentNode->id());
		if (!_subManager->isEmpty(currentNode->id())) {
			_subManager->undoSubs(currentNode->id(), *_rootSequence, _rateVec, _sp.get());
		}
	}
}

void rateMatrixSim::mutateSeqAlongBranch(tree::nodeP currentNode, int seqLength) {
	// std::cout << "mutating sequence along branch\n";
	// const int nodeId = currentNode->id();
	// const int parentId = currentNode->father()->id();
	const MDOUBLE distToFather = currentNode->dis2father();

	// std::cout << parentId << "->" << nodeId << "\n";
	// std::cout << "current branch length: " << distToFather << "\n";
	// mutateEntireSeq(currentNode, seqLength);

	// std::cin.get();
	if (distToFather > 2) {
		mutateEntireSeq(currentNode, seqLength);
	} else {
		mutateSeqGillespie(currentNode, seqLength, distToFather);
	}
	// _subManager->printSubManager();
}


void rateMatrixSim::mutateEntireSeq(tree::nodeP currentNode, int seqLength) {

	// std::cout << "mutating entire sequence!\n";
	const int nodeId = currentNode->id();
	const int parentId = currentNode->father()->id();


	for (size_t site = 0; site < seqLength; ++site) {
		int parentChar = (*_subManager).getCharacter(parentId, site, *_rootSequence);
		int nextChar = giveRandomChar(parentChar, nodeId ,_rateVec[site]);
		if (nextChar != parentChar){
			_subManager->handleEvent(parentId, site, parentChar, _rateVec, _sp.get(), *_rootSequence);
			if (currentNode->isLeaf()) {
				_subManager->handleEvent(nodeId, site, nextChar, _rateVec, _sp.get(), *_rootSequence);
			}
			// (*_rootSequence)[site] = nextChar;
		}
		// std::cout << site << "\n";
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

		int mutatedSite = _siteSampler->drawSample() - 1;
		int parentChar = (*_subManager).getCharacter(parentId, mutatedSite, *_rootSequence);
		int nextChar = giveRandomChar(parentChar, nodeId ,_rateVec[mutatedSite]);
		_subManager->handleEvent(parentId, mutatedSite, parentChar, _rateVec, _sp.get(), *_rootSequence);
		if (currentNode->isLeaf()) {
			_subManager->handleEvent(nodeId, mutatedSite, nextChar, _rateVec, _sp.get(), *_rootSequence);
		}
		// should update substitution per site variable,
		lambdaParam = _subManager->getReactantsSum();

		if (waitingTime < 0) {
			std::cout << branchLength << " " << lambdaParam << " " << waitingTime << "\n";
			errorMsg::reportError("waiting time is negative :(");
		}
		branchLength = branchLength - waitingTime;

		std::exponential_distribution<double> distribution(lambdaParam);
		waitingTime = distribution(_mt_rand);

	}
}

	// SubsMap nodeToSubsVec;
	// nodeToBlockMap[currentNode->name()] = BlockTree(sequenceSize);
	// size_t nodePosition = 0;
	// while (!nodes.empty()) {
	// 	nodes.pop();
	// 	if (!currentNode->isLeaf()) {
	// 		for (auto node: currentNode->getSons()) {
	// 			nodes.push(node);
	// 		}
	// 	} else {
	// 		if (nodes.empty()) break;
	// 	}
	// 	currentNode = nodes.top();
	// 	sequenceSize = nodeToBlockMap[currentNode->father()->name()].length() - 1;
	// 	BlockTree blocks = simulateAlongBranch(sequenceSize, currentNode->dis2father(), nodePosition);
	// 	nodeToBlockMap[currentNode->name()] = blocks;

	// 	++nodePosition;
	// }
	
// 	_avgSubtitutionsPerSite = 0.0;
// 	for (int p=0 ; p < _et->getRoot()->getNumberOfSons() ; ++p) {
// 	  recursiveGenerateSpecificSeq(rateVec, seqLength, _et->getRoot()->getSon(p));
// 	}
// 	_avgSubtitutionsPerSite /= 1.0*seqLength;
// }


// void rateMatrixSim::generate_rates_continuous_gamma(const int seqLength,const MDOUBLE alpha, Vdouble rates)
// {
// 	rates.clear();
// 	rates.resize(seqLength);
// 	for (int h = 0; h < seqLength; h++)  {
// 	  rates[h] = talRandom::SampleGamma(alpha);
// 	}
// }

// void rateMatrixSim::generate_seq_continuous_gamma(int seqLength) {
// 	sequence justAseq(_alph);
// 	_simulatedSequences.resize(_et->getNodesNum(),justAseq);
// 	for (int i=0; i < _simulatedSequences.size(); ++i) {
// 		_simulatedSequences[i].resize(seqLength);
// 	}
// 	generateRootSeq(seqLength); 

// 	vector<MDOUBLE> rateVec(seqLength);
// 	MDOUBLE alpha= (static_cast<gammaDistribution*>(_sp->distr()))->getAlpha();
// 	for (int h = 0; h < seqLength; h++)  {
// 	  rateVec[h] = talRandom::SampleGamma(alpha);
// 	}
	
// 	_avgSubtitutionsPerSite = 0.0;
// 	for (int p=0 ; p < _et->getRoot()->getNumberOfSons() ; ++p) {
// 	  recursiveGenerateSpecificSeq(rateVec, seqLength, _et->getRoot()->getSon(p));
// 	}
// 	_avgSubtitutionsPerSite /= 1.0*seqLength;
// }
		
// void rateMatrixSim::generate_seqWithRateVectorNoStopCodon(const Vdouble& simRates, int seqLength)
// {
// 	if (_alph->size() != 4)
// 		errorMsg::reportError("generate_seqWithRateVectorNoStopCodon is applicable only for nucleotide process");
// 	if (seqLength %3 != 0)
// 		errorMsg::reportError("generate_seqWithRateVectorNoStopCodon: seqLenth should be a multiplicative of 3");
// 	if (simRates.size() != seqLength)
// 		errorMsg::reportError("generate_seqWithRateVectorNoStopCodon: the size of simRates should be identical to seqLenth");

// //	sequence justAseq(_alph);
// //	vector<sequence> simulatedSequences(_et.getNodesNum(),justAseq);
// 	vector<sequence> simulatedSequences;
// 	//generate three nucleotide positions at a time. Repeat each position if the generated sequences contain stop codon 
// 	Vdouble rateVec(3);
// 	bool bStopCodonFound = false;
// 	codon codonAlph;
// 	for (int p = 0; p < seqLength; p+=3)
// 	{
// 		rateVec[0] = simRates[p];
// 		rateVec[1] = simRates[p+1];
// 		rateVec[2] = simRates[p+2];
// 		//generate 3 nucleotide positions with no stop codon
// 		for (int loop = 0; loop < 1000; ++loop)
// 		{
// 			bStopCodonFound = false;
// 			generate_seqWithRateVector(rateVec, 3);
// 			for (int s = 0; s < _simulatedSequences.size(); ++s)
// 			{
// 				string codonStr = _simulatedSequences[s].toString();
// 				if (codonAlph.isStopCodon(codonStr))
// 				{
//                     bStopCodonFound = true;
// 					break;
// 				}
// 			}
// 			if (!bStopCodonFound)
// 				break;
// 		}
// 		if (bStopCodonFound)
// 			errorMsg::reportError("Could not generate a position without stop codon");
// 		//append positions to the positions generated so far
// 		if (p == 0)
// 			simulatedSequences = _simulatedSequences; //this will copy also the names of the sequences
// 		else
// 		{
//             for (int i = 0; i < simulatedSequences.size(); ++i)
//                 simulatedSequences[i] += _simulatedSequences[i];
// 		}
// 	}
// 	_simulatedSequences = simulatedSequences;
// }



// void rateMatrixSim::generate_seqWithRateVector(const Vdouble& rateVec, const int seqLength) {
// 	sequence justAseq(_alph);
// 	_simulatedSequences.resize(_et->getNodesNum(),justAseq);
// 	for (int i=0; i < _simulatedSequences.size(); ++i) {
// 		_simulatedSequences[i].resize(seqLength);
// 	}
// 	generateRootSeq(seqLength); 

// 	_avgSubtitutionsPerSite = 0.0;
// 	for (int p=0 ; p < _et->getRoot()->getNumberOfSons() ; ++p) {
// 	  recursiveGenerateSpecificSeq(rateVec,seqLength,_et->getRoot()->getSon(p));
// 	}
// 	_avgSubtitutionsPerSite /= 1.0*seqLength;
// }

// void rateMatrixSim::generateRootSeq(int seqLength) {	
// 	for (int i = 0; i < seqLength; i++) {
// 		_simulatedSequences[_et->getRoot()->id()][i] =  giveRandomChar();
//      }

// 	_simulatedSequences[_et->getRoot()->id()].setAlphabet(_alph);
// 	_simulatedSequences[_et->getRoot()->id()].setName(_et->getRoot()->name());
// 	_simulatedSequences[_et->getRoot()->id()].setID(_et->getRoot()->id());

// }

void rateMatrixSim::generateRootSeq(int seqLength) {	
	for (int i = 0; i < seqLength; i++) {
		ALPHACHAR newChar = giveRandomChar();
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

// void rateMatrixSim::generateRootLog(int seqLength) {
// 	_rootSequence.resize(seqLength);

// 	for (int i = 0; i < seqLength; i++) {
// 		_rootSequence[i] =  giveRandomChar();
//      }

// 	_rootSequence.setAlphabet(_alph);
// 	_rootSequence.setName(_et->getRoot()->name());
// 	_rootSequence.setID(_et->getRoot()->id());
// }


// void rateMatrixSim::recursiveGenerateSpecificSeq(
// 							const vector<MDOUBLE> &rateVec,
// 							const int seqLength,
// 							tree::nodeP myNode) {

// 	for (int y = 0; y < seqLength; y++) {
// 		// MDOUBLE lenFromFather=myNode->dis2father()*rateVec[y];
// 		int aaInFather = _simulatedSequences[myNode->father()->id()][y];
// 		int newChar = giveRandomChar(aaInFather, myNode->id(), rateVec[y]);
// 		if(newChar != aaInFather) _avgSubtitutionsPerSite += 1;
// 		_simulatedSequences[myNode->id()][y] = newChar;
//     }
// 	_simulatedSequences[myNode->id()].setAlphabet(_alph);
// 	_simulatedSequences[myNode->id()].setName(myNode->name());
// 	_simulatedSequences[myNode->id()].setID(myNode->id());
// 	for (int x =0 ; x < myNode->getNumberOfSons(); ++x) {
// 	  recursiveGenerateSpecificSeq(rateVec, seqLength, myNode->getSon(x));
// 	}
// }



ALPHACHAR rateMatrixSim::giveRandomChar() const {
	for (int loop =0 ;loop<100000 ;loop++) {

		MDOUBLE theRandNum = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
		MDOUBLE sum = 0.0;
		for (int j=0;j<_alphaSize;++j) {
			sum+=_sp->freq(j);
			if (theRandNum<sum) return j;
		}
	} 
	errorMsg::reportError("Could not give random character. The reason is probably that the P_i do not sum to one.");
	return 1;
}

ALPHACHAR rateMatrixSim::giveRandomChar(const ALPHACHAR letterInFatherNode,
								 const int nodeId,
								 const MDOUBLE rateCat) const {
	// assert(letterInFatherNode>=0);
	// assert(letterInFatherNode<_alphaSize);
	int randChar = _cpijGam.getRandomChar(rateCat, nodeId, letterInFatherNode);
	// std::cout << randChar << "\n";
	return randChar;
	// for (int loop =0 ;loop<100000 ;loop++) {
	// 	MDOUBLE theRandNum = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
	// 	MDOUBLE sum = 0.0;
	// 	for (int j=0;j<_alphaSize;++j) {
	// 		// sum+=_sp->Pij_t(letterInFatherNode,j, length);
	// 		sum+=_cpijGam.getPij(rateCat, nodeId, letterInFatherNode, j);
	// 		if (theRandNum<sum) return j;
	// 	}
	// }
	// errorMsg::reportError("Could not give random character. The reason is probably that the Pij_t do not sum to one.");
	// return 1;
}


int rateMatrixSim::getRandCategory(const int pos) const {
  MDOUBLE theRandNum = talRandom::giveRandomNumberBetweenZeroAndEntry(1);
  MDOUBLE sum = 0.0;
  for (int j=0;j<_sp->categories() ;++j) {
     sum+=_sp->ratesProb(j);
   if (theRandNum<sum) return j;
  }
  errorMsg::reportError(" error in function simulateTree::getRandCategory() ");// also quit the program
  return -1;
}

sequenceContainer rateMatrixSim::toSeqData() {
	sequenceContainer myseqData;
	for (int i=0; i < _simulatedSequences.size(); ++i) {
		myseqData.add(_simulatedSequences[i]);
	}
	return myseqData;
}

std::unique_ptr<sequenceContainer> rateMatrixSim::toSeqDataWithoutInternalNodes() {
	std::unique_ptr<sequenceContainer> myseqData = std::make_unique<sequenceContainer>();
	// sequenceContainer myseqData;
	for (int i=0; i < _simulatedSequences.size(); ++i) {
		tree::nodeP theCurNode = _et->findNodeByName(_simulatedSequences[i].name());
		if (theCurNode == NULL)
			errorMsg::reportError("could not find the specified name: " + _simulatedSequences[i].name());
		if (theCurNode->isInternal()) continue;
		myseqData->add(_simulatedSequences[i]);
	}
	return myseqData;
}
