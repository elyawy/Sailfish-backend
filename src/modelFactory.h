#ifndef ___MODEL_FACTORY
#define ___MODEL_FACTORY

#include <vector>
#include <memory>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "../libs/Phylolib/includes/stochasticProcess.h"
#include "../libs/Phylolib/includes/alphabet.h"
#include "../libs/Phylolib/includes/tree.h"
#include "../libs/Phylolib/includes/amino.h"
#include "../libs/Phylolib/includes/nucleotide.h"
#include "../libs/Phylolib/includes/codon.h"

#include "../libs/Phylolib/includes/eigenAccelerator.h"
#include "../libs/Phylolib/includes/trivialAccelerator.h"
#include "../libs/Phylolib/includes/gammaDistribution.h"
#include "../libs/Phylolib/includes/customDistribution.h"

#include "CategorySampler.h"
#include "allModels.h"
#include "nonRevAccelerator.h"

// wrapper for all the information about the substitution model:
// alphabet = aa/nc
// model = JC/GTR/WAG etc...
// matrix_parameters = a b c d e f
// gamma = alpha
// gamma_categories = 1,2,3,4...

enum alphabetCode {
    NULLCODE,
    NUCLEOTIDE,
    AMINOACID,
    CODON
};




class modelFactory {

public:
    modelFactory(): 
        _alphabet(alphabetCode::NUCLEOTIDE),
        _model(modelCode::NUCJC),  // Default model
        _alpha(1.0),
        _gammaCategories(1) {}

    void buildModel(alphabetCode alphabet, modelCode model,
                    const std::vector<MDOUBLE>& parameters = {},
                    const std::string& modelFilePath = "") {
        _alphabet = alphabet;
        _model = model;
        _parameters = parameters;
        _modelFilePath = modelFilePath;

        size_t alphabetSize = 0;
        if (alphabet == alphabetCode::NUCLEOTIDE) alphabetSize = 4;
        else if (alphabet == alphabetCode::AMINOACID) alphabetSize = 20;
        else if (alphabet == alphabetCode::CODON) alphabetSize = 61;

        bool isReversible = (model != modelCode::NONREVERSIBLE);

        std::vector<MDOUBLE> qMatrix(alphabetSize * alphabetSize, 0.0);
        std::vector<MDOUBLE> frequencies(alphabetSize, 0.0);

        switch (_model) {
            case modelCode::NUCJC:
                _cachedRepModel = std::make_unique<nucJC>();
                break;
            case modelCode::AAJC:
                _cachedRepModel = std::make_unique<aaJC>();
                break;
            case modelCode::GTR: {
                Vdouble frequencies_(_parameters.begin(), _parameters.begin() + 4);
                _cachedRepModel = std::make_unique<gtrModel>(
                    frequencies_, _parameters[4], _parameters[5], _parameters[6],
                    _parameters[7], _parameters[8], _parameters[9]);
                break;
            }
            case modelCode::HKY: {
                Vdouble inProbabilities(_parameters.begin(), _parameters.begin() + 4);
                _cachedRepModel = std::make_unique<hky>(inProbabilities, _parameters[4]);
                break;
            }
            case modelCode::TAMURA92:
                _cachedRepModel = std::make_unique<tamura92>(_parameters[0], _parameters[1]);
                break;
            case modelCode::WYANGMODEL: {
                codon* coAlpha = dynamic_cast<codon*>(getAlphabet());
                _cachedRepModel = std::make_unique<wYangModel>(_parameters[0], _parameters[1], true, coAlpha);
                break;
            }
            case modelCode::REVERSIBLE: {
                std::ifstream in(_modelFilePath);
                if (!in.is_open()) throw std::runtime_error("Could not open file");
                std::stringstream contents;
                char buffer;
                while (in.get(buffer)) {
                    if (buffer == '\"' || buffer == '\n') continue;
                    contents << buffer;
                }
                datMatrixString aminoFileString(contents.str().c_str());
                _cachedRepModel = std::make_unique<pupAll>(aminoFileString);
                break;
            }
            case modelCode::NONREVERSIBLE: {
                // no repModel — accelerator built directly from raw Q + frequencies below
                _cachedRepModel = nullptr;
                // Note: For NONREV parameters should contain the flattened NxN Q matrix followed by N frequencies
                if (_parameters.size() != alphabetSize*alphabetSize + alphabetSize) 
                    throw std::runtime_error("Incorrect number of parameters for NONREV model");
                qMatrix.assign(_parameters.begin(), _parameters.end() - alphabetSize);
                frequencies.assign(_parameters.end() - alphabetSize, _parameters.end());
                break;
            }
            default: {
                auto it = modelToDatMatrixHolder.find(_model);
                if (it == modelToDatMatrixHolder.end())
                    throw std::runtime_error("Unknown model code.");
                _cachedRepModel = std::make_unique<pupAll>(it->second, alphabetSize);
                break;
            }
        }   
        
        _cachedPij = acceleratorPicker(alphabetSize, isReversible, 
                                       _cachedRepModel.get(), 
                                       qMatrix, frequencies);
    }

    std::unique_ptr<pijAccelerator> acceleratorPicker(size_t alphabetSize, bool isReversible,
        replacementModel* repModel = nullptr,
        const std::vector<MDOUBLE>& qMatrix = {},
        const std::vector<MDOUBLE>& frequencies = {}) {


        if (!isReversible) {
            switch (alphabetSize) {
                case 4:
                    return std::make_unique<nonRevAccelerator<4>>(qMatrix, frequencies);
                case 20:
                    return std::make_unique<nonRevAccelerator<20>>(qMatrix, frequencies);
                case 61:
                    return std::make_unique<nonRevAccelerator<61>>(qMatrix, frequencies);
                default:
                    throw std::runtime_error("Unsupported alphabet size for non-reversible model");
            }
        } 

        if (alphabetSize == 4) {
            return std::make_unique<trivialAccelerator>(repModel);
        } else if (alphabetSize == 20) {
            return std::make_unique<eigenAccelerator<20>>(repModel);
        } else if (alphabetSize == 61) {
            return std::make_unique<eigenAccelerator<61>>(repModel);
        } else {
            throw std::runtime_error("Unsupported alphabet size for reversible model");
        }
        
    }

    void setSiteRateModel(const std::vector<MDOUBLE>& rates,
                         const std::vector<MDOUBLE>& stationaryProbs,
                         const std::vector<std::vector<MDOUBLE>>& transitionMatrix = {}) {
        _customRates = rates;
        _rateCategoryProbs = stationaryProbs;
        _transitionMatrix = transitionMatrix;
    }

    std::vector<std::vector<MDOUBLE>> getEffectiveTransitionMatrix() const {
        if (_transitionMatrix.empty()) {
            // Independent rates: P[i][j] = π[j]
            size_t n = _rateCategoryProbs.size();
            return std::vector<std::vector<MDOUBLE>>(n, _rateCategoryProbs);
        }
        return _transitionMatrix;
    }

    const std::vector<MDOUBLE>& getRateCategoryProbs() const {
        return _rateCategoryProbs;
    }

    alphabet* getAlphabet() {
        // Return existing alphabet if already created
        if (_alphPtr) {
            return _alphPtr.get();
        }
        
        // Check if alphabet was set
        if (_alphabet == alphabetCode::NULLCODE) {
            std::cout << "alphabet was not set! returning null pointer\n";
            return nullptr;
        }
        
        // Create alphabet based on type
        if (_alphabet == alphabetCode::NUCLEOTIDE) {
            _alphPtr = std::make_unique<nucleotide>();
        } else if (_alphabet == alphabetCode::AMINOACID) {
            _alphPtr = std::make_unique<amino>();
        } else if (_alphabet == alphabetCode::CODON)
        {
            _alphPtr = std::make_unique<codon>();
        } else {
            return nullptr;
        }
        
        return _alphPtr.get();
    }

    // Get stochastic process using cached replacement model and current rate model
    // This is cheap - can be called many times with different rate models
    std::shared_ptr<stochasticProcess> getStochasticProcess() {
        if (!_cachedPij) {
            throw std::runtime_error("Model not built — call buildModel() first.");
        }

        // Create distribution with current rate model parameters
        _cachedDist = std::make_unique<customDistribution>(_customRates, _rateCategoryProbs);

        return std::make_shared<stochasticProcess>(_cachedDist.get(), _cachedPij.get());
    }

    CategorySampler getRateCategorySampler(size_t maxPathLength = 0) {

        // Create category sampler with current rate model parameters and return
        auto transitionMatrix = getEffectiveTransitionMatrix();
        return CategorySampler(transitionMatrix, _rateCategoryProbs, maxPathLength);
    }

    ~modelFactory() {}

private:
    std::unique_ptr<alphabet> _alphPtr;
    alphabetCode _alphabet;
    modelCode _model;
    std::string _modelFilePath;
    std::vector<MDOUBLE> _parameters;
    MDOUBLE _alpha;  // Kept for potential future use, not currently used
    size_t _gammaCategories;  // Kept for potential future use, not currently used
    std::vector<MDOUBLE> _customRates;
    std::vector<std::vector<MDOUBLE>> _transitionMatrix;
    std::vector<MDOUBLE> _rateCategoryProbs; // (stationary) probabilities over rate categories — NOT nucleotide/amino acid frequencies

    // Cached replacement model components (expensive to build)
    std::unique_ptr<replacementModel> _cachedRepModel; // Defines the Q matrix and equilibrium frequencies
    std::unique_ptr<pijAccelerator> _cachedPij;
    std::unique_ptr<customDistribution> _cachedDist; // Defines the distribution over rate categories for the stochastic process
};


#endif