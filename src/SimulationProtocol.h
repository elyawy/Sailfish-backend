#ifndef SIMULATIONPROTOCOL_H
#define SIMULATIONPROTOCOL_H

#include <sstream>

#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/DiscreteDistribution.h"



class SimulationProtocol
{
private:
    tree* _tree;
    size_t _numberOfBranches;
    size_t _sequenceSize;
    size_t _minSequenceSize;
    size_t _seed;
    std::vector<DiscreteDistribution*> _insertionLengthDistributions;
    std::vector<DiscreteDistribution*> _deletionLengthDistributions;
    std::vector<double> _insertionRates;
    std::vector<double> _deletionRates;
    bool _isSaveAncestral;

public:
    SimulationProtocol(tree* phylotree) : _tree(phylotree),
                                          _numberOfBranches(phylotree->getNodesNum()-1),
                                          _minSequenceSize(0) {}

    void setInsertionLengthDistributions(std::vector<DiscreteDistribution*> lengthDistributions) {
        if (lengthDistributions.size() != _numberOfBranches) {
            std::stringstream errorstr;
            errorstr << "Number of insertion length distributions and branches mismatch:\n";
            errorstr << "lengthDistributions.size()=" << lengthDistributions.size() << "\n";
            errorstr << "_numberOfBranches=" << _numberOfBranches << "\n";
            errorMsg::reportError(errorstr.str());
        }
        _insertionLengthDistributions = lengthDistributions;
    }

    DiscreteDistribution* getInsertionDistribution(size_t position) {
        if (_insertionLengthDistributions[position] == nullptr) {
            errorMsg::reportError("Null insertion length distribution accessed\n");
        }
        return _insertionLengthDistributions[position];
    }

    void setDeletionLengthDistributions(std::vector<DiscreteDistribution*> lengthDistributions) {
        if (lengthDistributions.size() != _numberOfBranches) {
            std::stringstream errorstr;
            errorstr << "Number of deletion legnth distributions and branches mismatch:\n";
            errorstr << "lengthDistributions.size()=" << lengthDistributions.size() << "\n";
            errorstr << "_numberOfBranches=" << _numberOfBranches << "\n";
            errorMsg::reportError(errorstr.str());
        }
        _deletionLengthDistributions = lengthDistributions;
    }

    DiscreteDistribution* getDeletionDistribution(size_t position) {
        if (_deletionLengthDistributions[position] == nullptr) {
            errorMsg::reportError("Null deletion length distribution accessed\n");
        }
        return _deletionLengthDistributions[position];
    }

    void setInsertionRates(std::vector<double> rates)
     {
        if (rates.size() != _numberOfBranches) {
            std::stringstream errorstr;
            errorstr << "Number of insertion rates and branches mismatch:\n";
            errorstr << "rates.size()=" << rates.size() << "\n";
            errorstr << "_numberOfBranches=" << _numberOfBranches << "\n";
            errorMsg::reportError(errorstr.str());
        }
        _insertionRates = rates;
    }

    double getInsertionRate(size_t position) {
        if (position >= _insertionRates.size()) {
            errorMsg::reportError("Null insertion rate accessed\n");
        }
        return _insertionRates[position];
    }

    void setDeletionRates(std::vector<double> rates) {
        if (rates.size() != _numberOfBranches) {
            std::stringstream errorstr;
            errorstr << "Number of deletion rates and branches mismatch:\n";
            errorstr << "rates.size()=" << rates.size() << "\n";
            errorstr << "_numberOfBranches=" << _numberOfBranches << "\n";
            errorMsg::reportError(errorstr.str());
        }
        _deletionRates = rates;
    }

    double getDeletionRate(size_t position) {
        if (position >= _deletionRates.size()) {
            errorMsg::reportError("Null deletion rate accessed\n");
        }
        return _deletionRates[position];
    }

    void setSequenceSize(size_t sequenceSize) {
        _sequenceSize = sequenceSize;
    }

    void setMinSequenceSize(size_t minSequenceSize) {
        _minSequenceSize = minSequenceSize;
    }

    void setSeed(size_t seed) {
        // Golden Ratio constant used for better hash scattering
        // Makes close by seeds produce very different _seed values
        // which reduces correlation when used in multiple simulations with close by seeds
        uint64_t phi = 0x9e3779b97f4a7c15;
        _seed = seed*phi;
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

    size_t getMinSequenceSize(){
        return _minSequenceSize;
    }

    void setSaveAncestral(bool saveAncestral) {
        _isSaveAncestral = saveAncestral;
    }

    bool getSaveAncestral() {
        return _isSaveAncestral;
    }




    ~SimulationProtocol() {}
};

#endif // SIMULATIONPROTOCOL_H