#ifndef BRANCH_TRANSITION_PROBABILITIES_H
#define BRANCH_TRANSITION_PROBABILITIES_H

#include <vector>
#include "../libs/Phylolib/includes/DiscreteNDistribution.h"
#include "../libs/Phylolib/includes/stochasticProcess.h"


constexpr uint8_t differenceToPowerOfTwo(int alphabetSize) {
    return (1 << static_cast<uint8_t>(std::ceil(std::log2(alphabetSize)))) - alphabetSize;
}


template<size_t AlphabetSize>
class BranchTransitionProbabilities {
public:
    BranchTransitionProbabilities(const double branchLength, const stochasticProcess& _sp)
    {
        const size_t numCategories = _sp.categories();
        _distributions.reserve(numCategories * AlphabetSize);

        for (size_t cat = 0; cat < numCategories; ++cat) {
            const double rate = _sp.rates(cat);
            
            for (size_t i = 0; i < AlphabetSize; ++i) {
                std::vector<double> probabilities(AlphabetSize+differenceToPowerOfTwo(AlphabetSize), 0.0);
                
                for (size_t j = 0; j < AlphabetSize; ++j) {
                    double prob = _sp.Pij_t(i, j, branchLength * rate);

                    probabilities[j] = (prob);
                }

                
                _distributions.emplace_back((probabilities));
            }
        }
    }
    

    DiscreteNDistribution<AlphabetSize+differenceToPowerOfTwo(AlphabetSize)>& getDistribution(int category, int character)  {
        size_t distributionIndex = category * AlphabetSize + character;

        return _distributions[distributionIndex];
    }


private:
    std::vector<DiscreteNDistribution<AlphabetSize+differenceToPowerOfTwo(AlphabetSize)>> _distributions;

};

#endif // BRANCH_TRANSITION_PROBABILITIES_H