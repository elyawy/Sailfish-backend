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
    BranchTransitionProbabilities(const double branchLength, const stochasticProcess& _sp, bool useVoseAlias = true)
        : _useVoseAlias(useVoseAlias)
    {
        const size_t numCategories = _sp.categories();
        constexpr size_t PaddedSize = nearestHigherPowerOfTwo(AlphabetSize);
        std::array<double, PaddedSize> probabilities{};

        if (_useVoseAlias) {
            _distributions.reserve(numCategories * AlphabetSize);
            for (size_t cat = 0; cat < numCategories; ++cat) {
                const double rate = _sp.rates(cat);
                for (size_t i = 0; i < AlphabetSize; ++i) {
                    probabilities.fill(0.0);
                    for (size_t j = 0; j < AlphabetSize; ++j) {
                        probabilities[j] = _sp.Pij_t(i, j, branchLength * rate);
                    }
                    _distributions.emplace_back(probabilities);
                }
            }
        } else {
            _cdfs.reserve(numCategories * AlphabetSize);
            for (size_t cat = 0; cat < numCategories; ++cat) {
                const double rate = _sp.rates(cat);
                for (size_t i = 0; i < AlphabetSize; ++i) {
                    std::array<double, AlphabetSize> cdf{};
                    double cumsum = 0.0;
                    for (size_t j = 0; j < AlphabetSize; ++j) {
                        cumsum += _sp.Pij_t(i, j, branchLength * rate);
                        cdf[j] = cumsum;
                    }
                    _cdfs.push_back(cdf);
                }
            }
        }
    }

    const DiscreteNDistribution<nearestHigherPowerOfTwo(AlphabetSize)>& getDistribution(int category, int character) const {
        return _distributions[category * AlphabetSize + character];
    }

    template<typename RngType>
    int drawSample(int category, int character, RngType& rng) const {
        if (_useVoseAlias) {
            return getDistribution(category, character).drawSample(rng);
        } else {
            const auto& cdf = _cdfs[category * AlphabetSize + character];
            // use single rng call to generate a uniform random number in [0, 1)
            double u = static_cast<double>(rng()) * (1.0 / static_cast<double>(std::numeric_limits<uint64_t>::max()));
            for (size_t j = 0; j < AlphabetSize; ++j) {
                if (u < cdf[j]) return static_cast<int>(j) + 1;
            }
            return static_cast<int>(AlphabetSize); // fallback: last character
        }
    }

private:
    bool _useVoseAlias;
    std::vector<DiscreteNDistribution<nearestHigherPowerOfTwo(AlphabetSize)>> _distributions;
    std::vector<std::array<double, AlphabetSize>> _cdfs;
};

#endif // BRANCH_TRANSITION_PROBABILITIES_H