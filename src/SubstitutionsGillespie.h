#ifndef ___SUBSTITUTIONS_GILLESPIE
#define ___SUBSTITUTIONS_GILLESPIE

#include <vector>
#include <random>
#include <limits>
#include <memory>

#include "../libs/Phylolib/includes/definitions.h"
#include "../libs/Phylolib/includes/stochasticProcess.h"
#include "../libs/Phylolib/includes/sequence.h"
#include "../libs/Phylolib/includes/DiscreteNDistribution.h"
#include "FastRejectionSampler.h"

template<typename RngType, size_t AlphabetSize>
class SubstitutionGillespie {
private:
    const stochasticProcess* _sp;
    RngType* _rng;
    std::vector<std::unique_ptr<DiscreteNDistribution<AlphabetSize>>> _gillespieSampler;
    MDOUBLE _minWeight;
    MDOUBLE _maxWeight;
    std::unique_ptr<FastRejectionSampler> _siteSampler;

    void initializeCharacterSamplers() {
        _gillespieSampler.resize(AlphabetSize);
        for (size_t i = 0; i < AlphabetSize; ++i) {
            std::vector<double> qRates(AlphabetSize, 0.0);
            double sum = -_sp->Qij(i, i);
            double normalizer = 1.0 / sum;
            for (size_t j = 0; j < AlphabetSize; ++j) {
                if (i == j) continue;
                qRates[j] = _sp->Qij(i, j) * normalizer;
            }
            _gillespieSampler[i] = std::make_unique<DiscreteNDistribution<AlphabetSize>>(qRates);
        }
    }

    void computeWeightBounds() {
        MDOUBLE minQii = std::numeric_limits<MDOUBLE>::max();
        MDOUBLE maxQii = 0.0;
        
        for (size_t i = 0; i < AlphabetSize; ++i) {
            MDOUBLE qii = -_sp->Qij(i, i);
            minQii = std::min(minQii, qii);
            maxQii = std::max(maxQii, qii);
        }
        
        MDOUBLE minRate = std::numeric_limits<MDOUBLE>::max();
        MDOUBLE maxRate = 0.0;
        
        for (size_t i = 0; i < _sp->categories(); ++i) {
            MDOUBLE rate = _sp->rates(i);
            minRate = std::min(minRate, rate);
            maxRate = std::max(maxRate, rate);
        }
        
        _minWeight = (minRate * minQii) / 2.0;
        _maxWeight = (maxRate * maxQii) * 2.0;
    }

public:
    SubstitutionGillespie(const stochasticProcess* sp, RngType* rng)
        : _sp(sp), _rng(rng), _siteSampler(nullptr) {
        initializeCharacterSamplers();
        computeWeightBounds();
    }

    void initializeSiteSampler(const sequence& seq, const std::vector<size_t>& rateCategories) {
        const size_t seqLength = seq.seqLen();
        std::vector<double> siteWeights(seqLength);
        for (size_t site = 0; site < seqLength; ++site) {
            ALPHACHAR currentChar = seq[site];
            MDOUBLE rate = _sp->rates(rateCategories[site]);
            siteWeights[site] = -_sp->Qij(currentChar, currentChar) * rate;
        }
        if (!_siteSampler) {
            _siteSampler = std::make_unique<FastRejectionSampler>(siteWeights, _minWeight, _maxWeight);
            return;
        } else {
            std::cout << "updating site sampler weights\n";
            _siteSampler->updateWeightBulk(siteWeights);
        }
    }

    /**
     * Simulate substitutions along a branch using Gillespie algorithm
     * 
     * @param currentSequence - sequence to mutate (modified in place)
     * @param branchLength - evolutionary time
     * @param rateCategories - gamma rate category for each site
     */
    void mutate(
        sequence& currentSequence,
        MDOUBLE branchLength,
        const std::vector<size_t>& rateCategories,
        bool rebuildGillespieSiteSampler)
    {
        const size_t seqLength = currentSequence.seqLen();
                
        // Initialize fast rejection sampler for this branch if needed
        if (!_siteSampler || rebuildGillespieSiteSampler) {
            initializeSiteSampler(currentSequence, rateCategories);
        }

        // Gillespie loop
        MDOUBLE totalTime = 0.0;
        while (totalTime < branchLength) {
            // Sample waiting time
            MDOUBLE lambda = _siteSampler->getSumOfWeights();
            std::exponential_distribution<double> exp_dist(lambda);
            MDOUBLE dt = exp_dist(*_rng);
            totalTime += dt;
            
            if (totalTime >= branchLength) break;
            
            // Sample which site mutates (proportional to rate)
            size_t mutatedSite = _siteSampler->sample(*_rng);
            ALPHACHAR oldChar = currentSequence[mutatedSite];
            
            // Sample new character from Q matrix row
            ALPHACHAR newChar = _gillespieSampler[oldChar]->drawSample(*_rng) - 1;
            currentSequence[mutatedSite] = newChar;
            
            // Update weight for mutated site
            MDOUBLE newWeight = -_sp->Qij(newChar, newChar) * _sp->rates(rateCategories[mutatedSite]);
            _siteSampler->updateWeight(mutatedSite, newWeight);
        }
    }
};

#endif // ___SUBSTITUTIONS_GILLESPIE