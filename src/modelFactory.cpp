#include "../libs/Phylolib/includes/amino.h"
#include "../libs/Phylolib/includes/nucleotide.h"
#include "../libs/Phylolib/includes/chebyshevAccelerator.h"
#include "../libs/Phylolib/includes/gammaDistribution.h"


#include "allModels.h"
#include "modelFactory.h"




void modelFactory::setAlphabet(alphabetCode alphabet) {
    if (_state!=factoryState::ALPHABET) {
        std::cout << "Please reset model if you wish to change alphabet.\n";
        return;
    }
    _alphabet = alphabet;
    _state = factoryState::MODEL;
}


void modelFactory::setReplacementModel(modelCode model) {
    if (_state!=factoryState::MODEL) {
        std::cout << "Please specify an alphabet before selecting a model.\n";
        return;
    }
    _model = model;
    _state = factoryState::PARAMETERS;
    if (_alphabet == alphabetCode::AMINOACID) _state = factoryState::GAMMA;
    if (_model == modelCode::AAJC || _model == modelCode::NUCJC) _state = factoryState::GAMMA;

}

void modelFactory::setModelParameters(std::vector<MDOUBLE> params) {
    if (_state!=factoryState::PARAMETERS) {
        std::cout << "Please specify an appropriate model before setting parameters.\n";
        return;
    }
    _parameters = params;
    _state = factoryState::GAMMA;
}

void modelFactory::setGammaParameters(MDOUBLE alpha, size_t numCategories) {
    if (_state!=factoryState::GAMMA) {
        std::cout << "Please specify a model before proceeding.\n";
        return;
    }
    _alpha = alpha;
    _gammaCategories = numCategories;
    _state = factoryState::COMPLETE;
}

void modelFactory::resetFactory() {
    _state = factoryState::ALPHABET;
}




std::unique_ptr<stochasticProcess> modelFactory::getStochasticProcess() {
    if (_state!=factoryState::COMPLETE) {
        std::cout << "Please set all the required model parameters.\n";
        // return;
    }

    alphabet *selectedAlphabet;
    nucleotide nuc;
    amino ami;
    if (_alphabet==alphabetCode::NUCLEOTIDE) selectedAlphabet = &nuc;
    if (_alphabet==alphabetCode::AMINOACID) selectedAlphabet = &ami;
    

    std::unique_ptr<replacementModel> repModel;


    if (_model==modelCode::NUCJC) repModel = std::make_unique<nucJC>();
    if (_model==modelCode::AAJC) repModel = std::make_unique<aaJC>();
    // if (_model==modelCode::GTR) repModel = std::make_unique<gtrModel>(_parameters);
    // if (_model==modelCode::HKY) repModel = &hky();
    // if (_model==modelCode::TAMURA92) repModel = &tamura92();
    // if (_model==modelCode::WYANGMODEL) repModel = &wYangModel();
    if (_model==modelCode::CPREV45) repModel = std::make_unique<pupAll>(datMatrixHolder::cpREV45);
    if (_model==modelCode::DAYHOFF) repModel = std::make_unique<pupAll>(datMatrixHolder::dayhoff);
    if (_model==modelCode::JONES) repModel = std::make_unique<pupAll>(datMatrixHolder::jones);
    if (_model==modelCode::MTREV24) repModel = std::make_unique<pupAll>(datMatrixHolder::mtREV24);
    if (_model==modelCode::WAG) repModel = std::make_unique<pupAll>(datMatrixHolder::wag);
    if (_model==modelCode::HIVB) repModel = std::make_unique<pupAll>(datMatrixHolder::HIVb);
    if (_model==modelCode::HIVW) repModel = std::make_unique<pupAll>(datMatrixHolder::HIVw);
    if (_model==modelCode::LG) repModel = std::make_unique<pupAll>(datMatrixHolder::lg);
    if (_model==modelCode::EMPIRICODON) repModel = std::make_unique<pupAll>(datMatrixHolder::empiriCodon);
    if (_model==modelCode::EX_BURIED) repModel = std::make_unique<pupAll>(datMatrixHolder::EX_BURIED);
    if (_model==modelCode::EX_EXPOSED) repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EXPOSED);
    if (_model==modelCode::EHO_EXTENDED) repModel = std::make_unique<pupAll>(datMatrixHolder::EHO_EXTENDED);
    if (_model==modelCode::EHO_HELIX) repModel = std::make_unique<pupAll>(datMatrixHolder::EHO_HELIX);
    if (_model==modelCode::EX_EHO_BUR_EXT) repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_BUR_EXT);
    if (_model==modelCode::EX_EHO_BUR_HEL) repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_EXP_HEL);
    if (_model==modelCode::EX_EHO_BUR_OTH) repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_BUR_OTH);
    if (_model==modelCode::EX_EHO_EXP_EXT) repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_EXP_EXT);
    if (_model==modelCode::EX_EHO_EXP_HEL) repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_EXP_HEL);
    if (_model==modelCode::EX_EHO_EXP_OTH) repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_EXP_OTH);

    chebyshevAccelerator accelerator(repModel.get());
    pijAccelerator* pij = &accelerator;

    gammaDistribution dist(_alpha, _gammaCategories);


    // stochasticProcess sp(&dist, pij);

    return std::make_unique<stochasticProcess>(&dist, pij);
}

alphabet* modelFactory::getAlphabet() {
    if (_alphabet==alphabetCode::NULLCODE) {
        std::cout << "alphabet was not set! returning null pointer\n";
        _alph = nullptr;
    }
    if (_alphabet==alphabetCode::NUCLEOTIDE) _alph = new nucleotide();
    if (_alphabet==alphabetCode::AMINOACID) _alph = new amino();
    return _alph;
}



modelFactory::~modelFactory()
{
    delete _alph;
}


