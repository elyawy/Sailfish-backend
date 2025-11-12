#ifndef ___CATEGORY_SAMPLER
#define ___CATEGORY_SAMPLER

#include <vector>
#include "../libs/Phylolib/includes/definitions.h"
#include "../libs/Phylolib/includes/errorMsg.h"
#include "../libs/Phylolib/includes/DiscreteDistribution.h"

/**
 * CategorySampler handles sampling rate categories with Markov chain autocorrelation.
 * 
 * Uses a transition matrix P[i][j] = probability of transitioning from category i to j.
 * The stationary distribution is used only for sampling the initial category.
 * 
 * Compatible with Yang (1995) auto-discrete-gamma model and extensions like G+I with autocorrelation.
 */
class CategorySampler {
public:
    /**
     * Constructor
     * @param transitionMatrix - Transition probability matrix P[i][j] (must be square, rows sum to 1)
     * @param stationaryProbs - Stationary distribution π for initial sampling (must sum to 1)
     */
    CategorySampler(const std::vector<std::vector<MDOUBLE>>& transitionMatrix,
                   const std::vector<MDOUBLE>& stationaryProbs)
        : _stationaryProbs(stationaryProbs), _previousCategory(-1) {
        
        // Validate inputs
        if (stationaryProbs.empty()) {
            errorMsg::reportError("CategorySampler: stationary probabilities cannot be empty");
        }
        
        if (transitionMatrix.empty()) {
            errorMsg::reportError("CategorySampler: transition matrix cannot be empty");
        }
        
        size_t numCategories = stationaryProbs.size();
        
        // Check matrix is square and matches stationary probs dimensions
        if (transitionMatrix.size() != numCategories) {
            errorMsg::reportError("CategorySampler: transition matrix dimensions don't match stationary probabilities");
        }
        
        MDOUBLE stationarySum = 0.0;
        for (size_t i = 0; i < numCategories; ++i) {
            if (transitionMatrix[i].size() != numCategories) {
                errorMsg::reportError("CategorySampler: transition matrix must be square");
            }
            
            // Validate row sums to 1 and non-negative
            MDOUBLE rowSum = 0.0;
            for (size_t j = 0; j < numCategories; ++j) {
                if (transitionMatrix[i][j] < 0.0) {
                    errorMsg::reportError("CategorySampler: transition probabilities must be non-negative");
                }
                rowSum += transitionMatrix[i][j];
            }
            
            if (std::abs(rowSum - 1.0) > 1e-6) {
                errorMsg::reportError("CategorySampler: each row of transition matrix must sum to 1");
            }
            
            // Validate stationary probs
            if (stationaryProbs[i] < 0.0) {
                errorMsg::reportError("CategorySampler: stationary probabilities must be non-negative");
            }
            stationarySum += stationaryProbs[i];
        }
        
        if (std::abs(stationarySum - 1.0) > 1e-6) {
            errorMsg::reportError("CategorySampler: stationary probabilities must sum to 1");
        }
        
        // Build transition samplers from provided matrix
        buildTransitionSamplers(transitionMatrix);
    }
    
    /**
     * Sample the next category
     * @return Category index
     */
    template<typename RngType = std::mt19937_64>
    int drawSample(RngType &rng) {
        if (_previousCategory < 0) {
            DiscreteDistribution initialSampler(_stationaryProbs);
            _previousCategory = initialSampler.drawSample(rng) - 1;
            return _previousCategory;
        }
        
        int nextCategory = _transitionSamplers[_previousCategory].drawSample(rng) - 1;
        _previousCategory = nextCategory;
        return nextCategory;
    }
    
    /**
     * Reset to sample from stationary distribution (for new sequences)
     */
    
    void reset() {
        // Sample new initial category from stationary distribution
        _previousCategory = -1;
    }

private:
    void buildTransitionSamplers(const std::vector<std::vector<MDOUBLE>>& transitionMatrix) {
        size_t numCategories = transitionMatrix.size();
        _transitionSamplers.clear();
        _transitionSamplers.reserve(numCategories);
        
        for (size_t i = 0; i < numCategories; ++i) {
            _transitionSamplers.emplace_back(transitionMatrix[i]);
        }
    }
    
    std::vector<MDOUBLE> _stationaryProbs;
    int _previousCategory;
    std::vector<DiscreteDistribution> _transitionSamplers;
};

#endif