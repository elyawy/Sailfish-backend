

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

    DiscreteDistribution* getInsertionDistribution(size_t position) const{
        return _insertionLengthDistributions[position];
    }

    std::vector<DiscreteDistribution*> getInsertionDistributions() const{
        return _insertionLengthDistributions;
    }

    void setDeletionLengthDistributions(std::vector<DiscreteDistribution*> lengthDistributions) {
        _deletionLengthDistributions = lengthDistributions;
    }

    DiscreteDistribution* getDeletionDistribution(size_t position) const{
        return _deletionLengthDistributions[position];
    }

    std::vector<DiscreteDistribution*> getDeletionDistributions() const{
        return _deletionLengthDistributions;
    }

    void setInsertionRates(std::vector<double> rates) {
        _insertionRates = rates;
    }

    double getInsertionRate(size_t position) const{
        return _insertionRates[position];
    }

    std::vector<double> getInsertionRates() const{
        return _insertionRates;
    }

    void setDeletionRates(std::vector<double> rates) {
        _deletionRates = rates;
    }

    double getDeletionRate(size_t position) const {
        return _deletionRates[position];
    }

    std::vector<double> getDeletionRates() const {
        return _deletionRates;
    }

    void setSequenceSize(size_t sequenceSize) {
        _sequenceSize = sequenceSize;
    }

    void setSeed(size_t seed) {
        _seed = seed;
    }

    size_t getSeed() const {
        return _seed;
    }

    tree* getTree() const {
        return _tree;
    }

    size_t getSequenceSize() const{
        return _sequenceSize;
    }

    void setSaveAncestral(bool saveAncestral) {
        _isSaveAncestral = saveAncestral;
    }

    bool getSaveAncestral() const {
        return _isSaveAncestral;
    }




    ~SimulationProtocol() {}
};
