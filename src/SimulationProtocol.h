#include <sstream>

#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/DiscreteDistribution.h"


enum class SiteRateModel {
    SIMPLE,         // Fast, no rate categories per site tracked
    INDEL_AWARE      // Slower, tracks rate categories per site affected by indel events
};

class SimulationProtocol
{
private:
    size_t _numberOfBranches;
    size_t _sequenceSize;
    size_t _minSequenceSize;
    std::vector<DiscreteDistribution*> _insertionLengthDistributions;
    std::vector<DiscreteDistribution*> _deletionLengthDistributions;
    std::vector<double> _insertionRates;
    std::vector<double> _deletionRates;
    SiteRateModel _siteRateModel;
    size_t _maxInsertionLength;

public:
    SimulationProtocol(size_t numberOfBranches) : 
    _numberOfBranches(numberOfBranches),
    _sequenceSize(0), _minSequenceSize(0),
    _siteRateModel(SiteRateModel::SIMPLE),
    _maxInsertionLength(0) {}

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

    void setGlobalInsertionLengthDistribution(DiscreteDistribution* lengthDistribution) {
        _insertionLengthDistributions.resize(_numberOfBranches, lengthDistribution);
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

    void setGlobalDeletionLengthDistribution(DiscreteDistribution* lengthDistribution) {
        _deletionLengthDistributions.resize(_numberOfBranches, lengthDistribution);
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

    void setGlobalInsertionRate(double rate) {
        _insertionRates.resize(_numberOfBranches, rate);
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

    void setGlobalDeletionRate(double rate) {
        _deletionRates.resize(_numberOfBranches, rate);
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

    size_t getSequenceSize(){
        return _sequenceSize;
    }
    
    void setMinSequenceSize(size_t minSequenceSize) {
        _minSequenceSize = minSequenceSize;
    }

    size_t getMinSequenceSize(){
        return _minSequenceSize;
    }


    void setIndelRateModel(SiteRateModel model) {
        _siteRateModel = model;
    }

    SiteRateModel getSiteRateModel() const {
        return _siteRateModel
    }

    void setMaxInsertionLength(size_t len) {
         _maxInsertionLength = len;
    }

    size_t getMaxInsertionLength() const {
        return _maxInsertionLength;
    }


    ~SimulationProtocol() {}
};
