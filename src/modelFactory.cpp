#include "../libs/Phylolib/includes/amino.h"
#include "../libs/Phylolib/includes/nucleotide.h"
#include "../libs/Phylolib/includes/chebyshevAccelerator.h"
#include "../libs/Phylolib/includes/trivialAccelerator.h"
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
    if (_model == modelCode::CUSTOM) _state = factoryState::MODEL_FILE;
}

void modelFactory::setModelParameters(std::vector<MDOUBLE> params) {
    if (_state!=factoryState::PARAMETERS) {
        std::cout << "Please specify an appropriate model before setting parameters.\n";
        return;
    }

    switch (_model){
        case modelCode::GTR: {
            if (params.size() != 10) {std::cout << "The 'GTR' model requires 10 paraemeters, " 
                << params.size() <<  " were provided\n"; return;}
            break;
        }    
        case modelCode::HKY: {
            if (params.size() != 5) {std::cout << "The 'HKY' model requires 5 paraemeters, " 
                << params.size() <<  " were provided\n"; return;}
            break;
        }    
        case modelCode::TAMURA92: {
            if (params.size() != 2) {std::cout << "The 'TAMURA92' model requires 2 paraemeters, " 
                << params.size() <<  " were provided\n"; return;}
            break;
        }    
        default:
            break;
    }

    _parameters = params;
    _state = factoryState::GAMMA;
}

void modelFactory::setCustomAAModelFile(const std::string &fileName) {
    if (_state!=factoryState::MODEL_FILE) {
        std::cout << "Please set the model to 'CUSTOM' before proceeding.\n";
        return;
    }
    _modelFilePath = fileName;
    _state = factoryState::GAMMA;
}

void modelFactory::setGammaParameters(MDOUBLE alpha, size_t numCategories) {
    if (_state!=factoryState::GAMMA) {
        std::cout << "Please specify a model and its correct parameters before proceeding.\n";
        return;
    }
    _alpha = alpha;
    _gammaCategories = numCategories;
    _state = factoryState::COMPLETE;
}

void modelFactory::resetFactory() {
    _state = factoryState::ALPHABET;
    _invariantProportion = 0.0;
}


std::shared_ptr<stochasticProcess> modelFactory::getStochasticProcess() {
    if (_state!=factoryState::COMPLETE) {
        std::cout << "Please set all the required model parameters.\n";
        // return;
    }    

    std::unique_ptr<replacementModel> repModel;


    switch (_model) {
        case modelCode::NUCJC:
            repModel = std::make_unique<nucJC>();
            break;
        case modelCode::AAJC:
            repModel = std::make_unique<aaJC>();
            break;
        case modelCode::GTR: {
            Vdouble frequencies = std::vector<MDOUBLE>(_parameters.begin(), _parameters.begin()+4);
            const MDOUBLE a2c = _parameters[4];
            const MDOUBLE a2g = _parameters[5];
            const MDOUBLE a2t = _parameters[6];
            const MDOUBLE c2g = _parameters[7];
            const MDOUBLE c2t = _parameters[8];
            const MDOUBLE g2t = _parameters[9];
            repModel = std::make_unique<gtrModel>(frequencies, a2c, a2g, a2t, c2g, c2t, g2t);
            break;
        }
        case modelCode::HKY: {
            Vdouble inProbabilities = std::vector<MDOUBLE>(_parameters.begin(), _parameters.begin()+4);
            const MDOUBLE TrTv = _parameters[4];
            repModel = std::make_unique<hky>(inProbabilities, TrTv); // handle different parameters!
            break;
        }
        case modelCode::TAMURA92: {
            const MDOUBLE theta = _parameters[0];
            const MDOUBLE TrTv = _parameters[1];
            repModel = std::make_unique<tamura92>(theta, TrTv);
            break;
        }
        case modelCode::WYANGMODEL:
            throw std::runtime_error("Model not implemented: " + std::to_string(static_cast<int>(_model)));
            break;
        case modelCode::CPREV45:
            repModel = std::make_unique<pupAll>(datMatrixHolder::cpREV45);
            break;
        case modelCode::DAYHOFF:
            repModel = std::make_unique<pupAll>(datMatrixHolder::dayhoff);
            break;
        case modelCode::JONES:
            repModel = std::make_unique<pupAll>(datMatrixHolder::jones);
            break;
        case modelCode::MTREV24:
            repModel = std::make_unique<pupAll>(datMatrixHolder::mtREV24);
            break;
        case modelCode::WAG:
            repModel = std::make_unique<pupAll>(datMatrixHolder::wag);
            break;
        case modelCode::HIVB:
            repModel = std::make_unique<pupAll>(datMatrixHolder::HIVb);
            break;
        case modelCode::HIVW:
            repModel = std::make_unique<pupAll>(datMatrixHolder::HIVw);
            break;
        case modelCode::LG:
            repModel = std::make_unique<pupAll>(datMatrixHolder::lg);
            break;
        case modelCode::EMPIRICODON:
            repModel = std::make_unique<pupAll>(datMatrixHolder::empiriCodon);
            break;
        case modelCode::EX_BURIED:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EX_BURIED);
            break;
        case modelCode::EX_EXPOSED:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EXPOSED);
            break;
        case modelCode::EHO_EXTENDED:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EHO_EXTENDED);
            break;
        case modelCode::EHO_HELIX:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EHO_HELIX);
            break;
        case modelCode::EHO_OTHER:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EHO_OTHER);
            break;
        case modelCode::EX_EHO_BUR_EXT:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_BUR_EXT);
            break;
        case modelCode::EX_EHO_BUR_HEL:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_EXP_HEL);
            break;
        case modelCode::EX_EHO_BUR_OTH:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_BUR_OTH);
            break;
        case modelCode::EX_EHO_EXP_EXT:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_EXP_EXT);
            break;
        case modelCode::EX_EHO_EXP_HEL:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_EXP_HEL);
            break;
        case modelCode::EX_EHO_EXP_OTH:
            repModel = std::make_unique<pupAll>(datMatrixHolder::EX_EHO_EXP_OTH);
            break;
        case modelCode::CUSTOM:
            std::ifstream in(_modelFilePath);
            if (!in.is_open()) throw std::runtime_error("Could not open file");
            std::stringstream contents;
            char buffer;
            while (in.get(buffer)) {
                if (buffer == '\"' || buffer == '\n') continue;
                contents << buffer;
            }
            in.close();
            const std::string &tmpstr = contents.str();
            const char* cstr = tmpstr.c_str();
            datMatrixString aminoFileString(cstr);
            repModel = std::make_unique<pupAll>(aminoFileString);
            break;
    }

    std::unique_ptr<pijAccelerator> pij;

    if (_alphabet==alphabetCode::AMINOACID) {
        pij = std::make_unique<chebyshevAccelerator>(repModel.get());
    } else if  (_alphabet==alphabetCode::NUCLEOTIDE){
        pij = std::make_unique<trivialAccelerator>(repModel.get());
    }

    gammaDistribution dist(_alpha, _gammaCategories);

    return std::make_shared<stochasticProcess>(&dist, pij.get());
}

alphabet* modelFactory::getAlphabet() {
    if (_alphabet==alphabetCode::NULLCODE) {
        std::cout << "alphabet was not set! returning null pointer\n";
        _alphPtr = nullptr;
    }
    // if (_alphabet==alphabetCode::NUCLEOTIDE) _alph = new nucleotide();
    // if (_alphabet==alphabetCode::AMINOACID) _alph = new amino();
    if (!_alphPtr) {
        if (_alphabet == alphabetCode::NUCLEOTIDE) {
            _alphPtr = std::make_unique<nucleotide>();
        } else if (_alphabet == alphabetCode::AMINOACID) {
            _alphPtr = std::make_unique<amino>();
        } else {
            return nullptr;
        }
    }
    return _alphPtr.get();
    // return _alph;
}

bool modelFactory::isModelValid() {
    return (_state==factoryState::COMPLETE);
}



modelFactory::~modelFactory() {}


