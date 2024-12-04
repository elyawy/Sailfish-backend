
#ifndef ___MODEL_FACTORY
#define ___MODEL_FACTORY

#include <vector>
#include <memory>
#include <sstream>

#include "../libs/Phylolib/includes/stochasticProcess.h"
#include "../libs/Phylolib/includes/alphabet.h"
#include "../libs/Phylolib/includes/tree.h"

// wrapper for all the information about the substitution model:
// alphabet = aa/nc
// model = JC/GTR/WAG etc...
// matrix_parameters = a b c d e f
// gamma = alpha
// gamma_categories = 1,2,3,4...

enum factoryState {
    ALPHABET,
    MODEL,
    PARAMETERS,
    MODEL_FILE,
    GAMMA,
    COMPLETE
};

enum alphabetCode {
    NULLCODE,
    NUCLEOTIDE,
    AMINOACID
};


enum modelCode {
    // nc:
    NUCJC,
    AAJC,
    GTR,
    HKY,
    TAMURA92,
    WYANGMODEL,
    // AA:
    CPREV45,
    DAYHOFF,
    JONES,	// THIS IS JTT
    MTREV24,
    WAG,
    HIVB,
    HIVW,
    LG,
    EMPIRICODON,
    EX_BURIED, 
    EX_EXPOSED,
    EHO_EXTENDED,
    EHO_HELIX,
    EHO_OTHER,
    EX_EHO_BUR_EXT,
    EX_EHO_BUR_HEL,
    EX_EHO_BUR_OTH,
    EX_EHO_EXP_EXT,
    EX_EHO_EXP_HEL,
    EX_EHO_EXP_OTH,
    CUSTOM
};



class modelFactory
{

public:
    modelFactory(tree* tr): _state(factoryState::ALPHABET), _tree(tr){
        _invariantProportion = 0.0;
    }

    std::shared_ptr<stochasticProcess> getStochasticProcess();

    void setAlphabet(alphabetCode alphabet);
    void setReplacementModel(modelCode model);
    void setModelParameters(std::vector<MDOUBLE> params);
    void setCustomAAModelFile(const std::string &fileName);
    void setGammaParameters(MDOUBLE alpha, size_t numCategories);
    void setInvariantSitesProportion(MDOUBLE invariantProportion) {_invariantProportion = invariantProportion;};

    MDOUBLE getInvariantSitesProportion() {return _invariantProportion;}

    void resetFactory();

    alphabet* getAlphabet();
    
    tree* getTree(){ return _tree;};

    bool isModelValid();

    ~modelFactory();
private:
    factoryState _state;
    tree* _tree;
    alphabet* _alph;
    alphabetCode _alphabet;
    modelCode _model;
    std::string _modelFilePath;
    std::vector<MDOUBLE> _parameters;
    MDOUBLE _alpha;
    size_t _gammaCategories;
    MDOUBLE _invariantProportion;

};


#endif
