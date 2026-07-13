
#ifndef ___ALL_MODELS
#define ___ALL_MODELS

// nucleotide models
#include "../libs/Phylolib/includes/nucJC.h"
#include "../libs/Phylolib/includes/aaJC.h"
#include "../libs/Phylolib/includes/gtrModel.h"
#include "../libs/Phylolib/includes/hky.h"
#include "../libs/Phylolib/includes/tamura92.h"
#include "../libs/Phylolib/includes/wYangModel.h"

// amino acid models
#include "../libs/Phylolib/includes/readDatMatrix.h"

enum modelCode {
    // nc:
    NUCJC,
    GTR,
    HKY,
    TAMURA92,
    // AA:
    AAJC, // amino acid JC
    CPREV45,
    DAYHOFF,
    JONES,	// THIS IS JTT
    MTREV24,
    WAG,
    HIVB,
    HIVW,
    LG,
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
    // CODON:
    EMPIRICODON,
    WYANGMODEL,
    // custom models:
    REVERSIBLE, // REVERSIBLE is for reversible models, and defined through file.
    NONREVERSIBLE, // NONREVERSIBLE is for non-reversible models, and defined through Q matrix and frequencies flat vector.
};

// get the model with a switch case
inline const datMatrixString& getDatMatrixStringForModel(modelCode model) {
    switch (model) {
        case CPREV45: return datMatrixHolder::cpREV45;
        case DAYHOFF: return datMatrixHolder::dayhoff;
        case JONES: return datMatrixHolder::jones;
        case MTREV24: return datMatrixHolder::mtREV24;
        case WAG: return datMatrixHolder::wag;
        case HIVB: return datMatrixHolder::HIVb;
        case HIVW: return datMatrixHolder::HIVw;
        case LG: return datMatrixHolder::lg;
        case EMPIRICODON: return datMatrixHolder::empiriCodon;
        case EX_BURIED: return datMatrixHolder::EX_BURIED;
        case EX_EXPOSED: return datMatrixHolder::EX_EXPOSED;
        case EHO_EXTENDED: return datMatrixHolder::EHO_EXTENDED;
        case EHO_HELIX: return datMatrixHolder::EHO_HELIX;
        case EHO_OTHER: return datMatrixHolder::EHO_OTHER;
        case EX_EHO_BUR_EXT: return datMatrixHolder::EX_EHO_BUR_EXT;
        case EX_EHO_BUR_HEL: return datMatrixHolder::EX_EHO_BUR_HEL;
        case EX_EHO_BUR_OTH: return datMatrixHolder::EX_EHO_BUR_OTH;
        case EX_EHO_EXP_EXT: return datMatrixHolder::EX_EHO_EXP_EXT;
        case EX_EHO_EXP_HEL: return datMatrixHolder::EX_EHO_EXP_HEL;
        case EX_EHO_EXP_OTH: return datMatrixHolder::EX_EHO_EXP_OTH;
        default:
            throw std::runtime_error("Unknown model code.");
    }
}


#endif