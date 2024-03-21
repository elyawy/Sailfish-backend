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

//class sequenceData; // to be able to go to simulate data.

class rateMatrixSim {
public:
	// explicit simulateTree(tree*  _inEt,const stochasticProcess* sp,
	// 	const alphabet* alph);
	explicit rateMatrixSim(modelFactory& mFac);

	// void generate_seq(int seqLength);
	void generate_substitution_log(int seqLength);
	void mutateSeqRecuresively(tree::nodeP currentNode, int seqLength);

	void setSeed(size_t seed);
	// int getSeed();

	// This function generates the sequences not using the discrete gamma, but rather,
	// the rates are sampled from the continuous distribution.
	// It assumes the Gamma distribution has mean 1 (alpha = beta).
	// void generate_seq_continuous_gamma(int seqLength);

	// void generate_seqWithRateVector(const Vdouble& simRates, const int seqLength);	
	//these function do the same simulation as above but check that no stop codon is created.
	//applicable only when the stochasticProcess is based on nucleotides
	// void generate_seqWithRateVectorNoStopCodon(const Vdouble& simRates, int seqLength);

	tree* gettree() {return _et;}
	virtual ~rateMatrixSim();
	sequenceContainer toSeqData();
	std::unique_ptr<sequenceContainer> toSeqDataWithoutInternalNodes();
	// void generate_rates_continuous_gamma(const int seqLength,const MDOUBLE alpha,Vdouble rates);
	MDOUBLE getAvgSub() {return _avgSubtitutionsPerSite;}
	
private:
	void generateRootSeq(int seqLength);

	// void recursiveGenerateSpecificSeq(const Vdouble& rateVec, int seqLength, tree::nodeP myNode);
	int giveRandomChar() const;
	int giveRandomChar(const int letterInFatherNode, const int nodeId,const MDOUBLE rateCat) const;
	int getRandCategory(const int pos) const;



	void mutateSeqAlongBranch(tree::nodeP parentNode, int seqLength);
	void mutateEntireSeq(tree::nodeP currentNode, int seqLength);
	void mutateSeqGillespie(tree::nodeP currentNode, int seqLength, MDOUBLE distToParent);
	void undoLastSubs(int fromNode);
	

	vector<sequence> _simulatedSequences; // the sequences (nodes * seqLen)
	std::shared_ptr<sequence> _rootSequence;
	vector<MDOUBLE> _rateVec;

	tree* _et;
	std::shared_ptr<const stochasticProcess> _sp;
	size_t _alphaSize;
	computePijGam _cpijGam;
	int _seed;
	const alphabet* _alph;
	MDOUBLE _avgSubtitutionsPerSite;
	MDOUBLE _subtitutionsRatePerSite;
	MDOUBLE _rateCategoriesSum;
	std::vector<MDOUBLE> _normalizedRateCategories;
	std::unique_ptr<DiscreteDistribution> _siteSampler;
	std::mt19937 _mt_rand;


	std::shared_ptr<substitutionManager> _subManager;
	
};

#endif

