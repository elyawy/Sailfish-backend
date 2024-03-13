

#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/DiscreteDistribution.h"



class SimulationProtocol
{
private:
    tree* _tree;
    size_t _sequenceSize;
    size_t _seed;
    std::vector<DiscreteDistribution*> _insertionLengthDistributions;
    std::vector<DiscreteDistribution*> _deletionLengthDistributions;
    std::vector<double> _insertionRates;
    std::vector<double> _deletionRates;
    bool _isSaveAncestral;

public:
    SimulationProtocol(tree* phylotree) : _tree(phylotree) {}

    void setInsertionLengthDistributions(std::vector<DiscreteDistribution*> lengthDistributions) {
        _insertionLengthDistributions = lengthDistributions;
    }

    DiscreteDistribution* getInsertionDistribution(size_t position) {
        return _insertionLengthDistributions[position];
    }

    void setDeletionLengthDistributions(std::vector<DiscreteDistribution*> lengthDistributions) {
        _deletionLengthDistributions = lengthDistributions;
    }

    DiscreteDistribution* getDeletionDistribution(size_t position) {
        return _deletionLengthDistributions[position];
    }

    void setInsertionRates(std::vector<double> rates) {
        _insertionRates = rates;
    }

    double getInsertionRate(size_t position) {
        return _insertionRates[position];
    }

    void setDeletionRates(std::vector<double> rates) {
        _deletionRates = rates;
    }

    double getDeletionRate(size_t position) {
        return _deletionRates[position];
    }

    void setSequenceSize(size_t sequenceSize) {
        _sequenceSize = sequenceSize;
    }

    void setSeed(size_t seed) {
        _seed = seed;
    }

    size_t getSeed() {
        return _seed;
    }

    tree* getTree(){
        return _tree;
    }

    size_t getSequenceSize(){
        return _sequenceSize;
    }

    void setSaveAncestral(bool saveAncestral) {
        _isSaveAncestral = saveAncestral;
    }

    bool getSaveAncestral() {
        return _isSaveAncestral;
    }




    ~SimulationProtocol() {}
};
