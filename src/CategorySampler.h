#ifndef ___CATEGORY_SAMPLER
#define ___CATEGORY_SAMPLER

#include <vector>
#include "../libs/Phylolib/includes/definitions.h"
#include "../libs/Phylolib/includes/errorMsg.h"
#include "../libs/Phylolib/includes/DiscreteDistribution.h"

/**
 * CategorySampler handles sampling rate categories with optional autocorrelation.
 * 
 * When correlation = 0: Categories are sampled independently from stationary distribution
 * When correlation > 0: Adjacent sites have correlated categories (autocorrelation model)
 * 
 * Uses transition matrix: P[i][j] = ρ * δ(i,j) + (1-ρ) * π[j]
 * where ρ is correlation and π is the stationary distribution
 */
class CategorySampler {
public:
    /**
     * Constructor
     * @param stationaryProbs - Stationary distribution π (probabilities of each category)
     * @param correlation - Autocorrelation parameter ρ ∈ [0, 1]
     */
    CategorySampler(const std::vector<MDOUBLE>& stationaryProbs, MDOUBLE correlation)
        : _stationaryProbs(stationaryProbs), _correlation(correlation), _previousCategory(-1) {
        
        if (correlation < 0.0 || correlation > 1.0) {
            errorMsg::reportError("CategorySampler: correlation must be between 0 and 1");
        }
        
        if (stationaryProbs.empty()) {
            errorMsg::reportError("CategorySampler: stationary probabilities cannot be empty");
        }
        
        buildTransitionMatrix();
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
    
    /**
     * Get current correlation parameter
     */
    MDOUBLE getCorrelation() const { return _correlation; }
    
    /**
     * Update correlation parameter and rebuild transition matrix
     */
    void setCorrelation(MDOUBLE correlation) {
        if (correlation < 0.0 || correlation > 1.0) {
            errorMsg::reportError("CategorySampler: correlation must be between 0 and 1");
        }
        
        if (_correlation == correlation) return;
        
        _correlation = correlation;
        buildTransitionMatrix();
        // Note: _previousCategory is not reset, allowing smooth transition
    }

private:
    void buildTransitionMatrix() {
        size_t numCategories = _stationaryProbs.size();
        _transitionSamplers.clear();
        _transitionSamplers.reserve(numCategories);
        
        for (size_t i = 0; i < numCategories; ++i) {
            std::vector<MDOUBLE> transitionProbs(numCategories);
            
            for (size_t j = 0; j < numCategories; ++j) {
                // P[i][j] = ρ * δ(i,j) + (1-ρ) * π[j]
                if (i == j) {
                    transitionProbs[j] = _correlation + (1.0 - _correlation) * _stationaryProbs[j];
                } else {
                    transitionProbs[j] = (1.0 - _correlation) * _stationaryProbs[j];
                }
            }
            
            _transitionSamplers.emplace_back(transitionProbs);
        }
    }
    
    MDOUBLE _correlation;
    std::vector<MDOUBLE> _stationaryProbs;
    int _previousCategory;
    std::vector<DiscreteDistribution> _transitionSamplers;
};

#endif