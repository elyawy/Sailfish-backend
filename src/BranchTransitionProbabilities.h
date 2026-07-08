#ifndef BRANCH_TRANSITION_PROBABILITIES_H
#define BRANCH_TRANSITION_PROBABILITIES_H

#include <vector>
#include "DiscreteNDistribution.h"
#include "../libs/Phylolib/includes/stochasticProcess.h"


constexpr uint8_t nearestHigherPowerOfTwo(int alphabetSize) {
    int p = 1;
    while (p < alphabetSize) p <<= 1;
    return static_cast<uint8_t>(p);
}


template<size_t AlphabetSize>
class BranchTransitionProbabilities {
public:
    BranchTransitionProbabilities(const double branchLength, const stochasticProcess& _sp)
    {
        const size_t numCategories = _sp.categories();
        _distributions.reserve(numCategories * AlphabetSize);
        constexpr size_t PaddedSize = nearestHigherPowerOfTwo(AlphabetSize);
        std::array<double, PaddedSize> probabilities{};

        for (size_t cat = 0; cat < numCategories; ++cat) {
            const double rate = _sp.rates(cat);
            
            for (size_t i = 0; i < AlphabetSize; ++i) {
                probabilities.fill(0.0);

                for (size_t j = 0; j < AlphabetSize; ++j) {
                    double prob = _sp.Pij_t(i, j, branchLength * rate);
                    probabilities[j] = (prob);
                }

                
                _distributions.emplace_back((probabilities));
            }
        }
    }
    

    const DiscreteNDistribution<nearestHigherPowerOfTwo(AlphabetSize)>& getDistribution(int category, int character) const {
        size_t distributionIndex = category * AlphabetSize + character;

        return _distributions[distributionIndex];
    }


private:
    std::vector<DiscreteNDistribution<nearestHigherPowerOfTwo(AlphabetSize)>> _distributions;

};

#endif // BRANCH_TRANSITION_PROBABILITIES_H