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

	// void setSeed(size_t seed);
	void setRng(mt19937_64 *rng);
	// const mt19937_64& getRng();
	void setNodesToSaves(std::vector<size_t> nodeIDs);
	void setSaveRates(bool saveRates);
	void clearRatesVec() { _siteRates.clear();}
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
	// sequenceContainer toSeqData();
	std::unique_ptr<sequenceContainer> getSequenceContainer();
	std::vector<double> getSiteRates() { return _siteRates;};
	// void generate_rates_continuous_gamma(const int seqLength,const MDOUBLE alpha,Vdouble rates);
	// MDOUBLE getAvgSub() {return _avgSubtitutionsPerSite;}
	
private:
	void generateRootSeq(int seqLength);
	void mutateSeqAlongBranch(tree::nodeP parentNode, int seqLength);
	void mutateEntireSeq(tree::nodeP currentNode, int seqLength);
	void mutateSeqGillespie(tree::nodeP currentNode, int seqLength, MDOUBLE distToParent);
	void saveSequence(const int &nodeId,const std::string &name);
	void initGillespieSampler();
	void setSaveStateLeaves(const tree::nodeP &node);
	bool testSumOfRates();
	// void resetSim();

	tree* _et;
	std::shared_ptr<const stochasticProcess> _sp;
	const alphabet* _alph;
	computePijGam _cpijGam;
	sequence _rootSequence;
	substitutionManager _subManager;
	std::vector<bool> _nodesToSave;
	bool _saveRates;
	std::uniform_real_distribution<double> _biased_coin;
	std::vector<std::unique_ptr<DiscreteDistribution>> _gillespieSampler;

	// vector<MDOUBLE> _rateVec;
	std::vector<size_t> _rateCategories;
	std::vector<double> _siteRates;
	std::unique_ptr<sequenceContainer> _simulatedSequences; // the sequences (nodes * seqLen)
	std::unique_ptr<DiscreteDistribution> _siteSampler;
	std::unique_ptr<DiscreteDistribution> _frequencySampler;
	std::unique_ptr<DiscreteDistribution> _rateSampler;

	std::mt19937_64 *_mt_rand;

	
};

#endif

