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
     * Sample a left-sided bridge: condition on left flanking category only
     * Used when insertion happens at the end of a block (no right neighbor)
     * @param leftCategory The rate category of the left flanking position
     * @param length Number of positions to sample
     * @return Vector of sampled rate categories
     */
    template<typename RngType = std::mt19937_64>
    std::vector<size_t> sampleLeftSidedBridge(size_t leftCategory, size_t length, RngType &rng) {
        // TODO: Implement proper bridge sampling conditioned on left category
        // For now, just sample forward from the left category
        std::vector<size_t> samples;
        samples.reserve(length);
        
        int previousCategory = static_cast<int>(leftCategory);
        for (size_t i = 0; i < length; ++i) {
            _previousCategory = previousCategory;
            int nextCategory = drawSample(rng);
            samples.push_back(nextCategory);
            previousCategory = nextCategory;
        }
        
        return samples;
    }

    /**
     * Sample a right-sided bridge: condition on right flanking category only
     * Used when insertion happens at the beginning of a block (no left neighbor)
     * @param rightCategory The rate category of the right flanking position
     * @param length Number of positions to sample
     * @return Vector of sampled rate categories
     */
    template<typename RngType = std::mt19937_64>
    std::vector<size_t> sampleRightSidedBridge(size_t rightCategory, size_t length, RngType &rng) {
        // TODO: Implement proper bridge sampling conditioned on right category
        // For now, just sample backward (reverse the process) from right category
        std::vector<size_t> samples;
        samples.reserve(length);
        
        // Sample forward and hope it ends near rightCategory (naive approach)
        // Proper implementation would use backward sampling
        _previousCategory = -1; // Start fresh
        for (size_t i = 0; i < length; ++i) {
            int nextCategory = drawSample(rng);
            samples.push_back(nextCategory);
        }
        
        return samples;
    }

    /**
     * Sample a bridge: condition on both left and right flanking categories
     * Used when insertion happens in the middle of a block
     * @param leftCategory The rate category of the left flanking position
     * @param rightCategory The rate category of the right flanking position
     * @param length Number of positions to sample
     * @return Vector of sampled rate categories
     */
    template<typename RngType = std::mt19937_64>
    std::vector<size_t> sampleBridge(size_t leftCategory, size_t rightCategory, size_t length, RngType &rng) {
        // TODO: Implement proper bridge sampling conditioned on both left and right
        // This requires a more sophisticated algorithm (e.g., forward-backward algorithm)
        // For now, just sample forward from left category
        std::vector<size_t> samples;
        samples.reserve(length);
        
        int previousCategory = static_cast<int>(leftCategory);
        for (size_t i = 0; i < length; ++i) {
            _previousCategory = previousCategory;
            int nextCategory = drawSample(rng);
            samples.push_back(nextCategory);
            previousCategory = nextCategory;
        }
        
        return samples;
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