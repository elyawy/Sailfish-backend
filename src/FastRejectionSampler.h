#include <vector>
#include <unordered_map>

#include <cmath>
#include <random>
#include <algorithm>

#include <iostream>

class FastRejectionSampler
{
private:

    std::vector<double> _weights;
    double _minWeight;
    double _maxWeight;
    std::uniform_real_distribution<double> _biasedCoin;
    double _totalWeightsSum;

    std::vector<std::vector<size_t>> _levelToWeights;
    std::unordered_map<size_t, size_t> _weightIndexToBin;

    int _minWeightLevel;
    int _maxWeightLevel;

    std::vector<double> _levelsWeights;

public:
    FastRejectionSampler(const std::vector<double> &weights, double minWeight, double maxWeight): 
        _weights(weights), _minWeight(minWeight), _maxWeight(maxWeight) , _biasedCoin(0.0,1.0), _totalWeightsSum(0.0) {

        _minWeightLevel = static_cast<int>(std::log2(minWeight));
        if (_minWeightLevel >= 0) _minWeightLevel += 1;
        _maxWeightLevel = static_cast<int>(std::log2(maxWeight)) + 1;
        size_t numLevels = _maxWeightLevel - _minWeightLevel + 1;

        _levelToWeights.resize(numLevels);
        _levelsWeights.resize(numLevels, 0.0);


        for(size_t i=0; i < _weights.size(); ++i) {
            MDOUBLE currentWeight = _weights[i];
            if (currentWeight == 0.0) continue;
            _totalWeightsSum += currentWeight;
            int level = static_cast<int>(std::log2(currentWeight));
            if (level >= 0) level += 1;
            level -= _minWeightLevel;
            _levelsWeights[level] += currentWeight;
            size_t innerIndex = _levelToWeights.at(level).size();
            _levelToWeights.at(level).push_back(i);
            _weightIndexToBin[i] = innerIndex;

        }



    }
    
    template <typename Generator>
    size_t sample(Generator&& gen) {
        double levelSampler = std::uniform_real_distribution<double>(0.0, _totalWeightsSum)(gen);
        int selectedLevel = 0;

        double cumulativeWeight = 0.0;

        for (size_t i = 0; i < _levelsWeights.size(); i++) {
            cumulativeWeight += _levelsWeights[i];
            selectedLevel = i;
            if (levelSampler < cumulativeWeight) break;
        }

        int correctedLevel =  selectedLevel + _minWeightLevel;
        double levelConversion = 1.0 / std::pow(2, correctedLevel);
        auto binsInSelectedLevel = _levelToWeights.at(selectedLevel);

        
        auto binSampler = std::uniform_int_distribution<int>(0, binsInSelectedLevel.size() - 1);

        // rejection sampling:
        while (true) {
            int selectedBin = binSampler(gen);
            size_t selectedIndex = binsInSelectedLevel[selectedBin];
            double rejection = _weights[selectedIndex] * levelConversion;
            double coinToss = _biasedCoin(gen);
            if (coinToss < rejection) return selectedIndex;
        }
    
        // return UINT64_MAX;
    }


    void updateWeight(int weightIndex, double newWeight) {
        if ((newWeight > _maxWeight) || (newWeight < _minWeight)) {
            std::cout << "new weight is out of bounds\n";
            abort();
        }
        double oldWeight = _weights[weightIndex];
        int oldLevel = static_cast<int>(std::log2(_weights[weightIndex]));
        size_t oldBinIndex = _weightIndexToBin.at(weightIndex);
        if (oldLevel >= 0) oldLevel += 1;
        int oldLevelIndex = oldLevel - _minWeightLevel;

        int newLevel = static_cast<int>(std::log2(newWeight));
        if (newLevel >= 0) newLevel += 1;
        int newLevelIndex = newLevel - _minWeightLevel;

        _totalWeightsSum -= oldWeight;
        _totalWeightsSum += newWeight;


        if (oldLevel == newLevel) {
            _levelsWeights[newLevelIndex] -= oldWeight;
            _levelsWeights[newLevelIndex] += newWeight;

            _weights[weightIndex] = newWeight;
            return;
        }

        // remove weight from old level
        _levelsWeights[oldLevelIndex] -= oldWeight;
        size_t numberOfBins = _levelToWeights.at(oldLevelIndex).size();
        if (numberOfBins > 1) {
            size_t weightIndexOfLastBin = _levelToWeights.at(oldLevelIndex)[numberOfBins - 1];
            _levelToWeights.at(oldLevelIndex)[oldBinIndex] = weightIndexOfLastBin;
            _weightIndexToBin[weightIndexOfLastBin] = oldBinIndex;
        } else {
            _levelsWeights[oldLevelIndex] = 0.0;
        }
        _levelToWeights.at(oldLevelIndex).pop_back();

        // add to new level
        _levelsWeights[newLevelIndex] += newWeight;
        _weightIndexToBin[weightIndex] = _levelToWeights.at(newLevelIndex).size();
        _levelToWeights.at(newLevelIndex).push_back(weightIndex);

        // update weight in the original weights vector
        _weights[weightIndex] = newWeight;

    }

    const std::vector<double> & getLevelsWeights() {
        return _levelsWeights;
    }


    double getLevelWeight(int level) {
        return _levelsWeights.at(level);
    }

    double getSumOfWeights() {
        return _totalWeightsSum;
    }

    int getLevelBin(int level) {
        return level + _minWeightLevel;
    }

    bool checkValidity(){
        double epsilon = 1e-10;
        double sum = 0.0;
        for(size_t i=0; i < _weights.size(); ++i) {
            sum += _weights[i];
        }
        if (abs(sum - _totalWeightsSum) > epsilon) return false;

        sum = 0.0;
        for (size_t i = 0; i < _levelsWeights.size(); i++) {
            sum += _levelsWeights[i];
        }
        if (abs(sum - _totalWeightsSum) > epsilon) return false;

        return true;
    }

    ~FastRejectionSampler(){};
};

