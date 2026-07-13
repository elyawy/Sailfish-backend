
#ifndef ___ALL_MODELS
#define ___ALL_MODELS


#include <unordered_map>

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

static const std::unordered_map<modelCode, datMatrixString> modelToDatMatrixHolder = {
    {modelCode::CPREV45, datMatrixHolder::cpREV45},
    {modelCode::DAYHOFF, datMatrixHolder::dayhoff},
    {modelCode::JONES, datMatrixHolder::jones},
    {modelCode::MTREV24, datMatrixHolder::mtREV24},
    {modelCode::WAG, datMatrixHolder::wag},
    {modelCode::HIVB, datMatrixHolder::HIVb},
    {modelCode::HIVW, datMatrixHolder::HIVw},
    {modelCode::LG, datMatrixHolder::lg},
    {modelCode::EMPIRICODON, datMatrixHolder::empiriCodon},
    {modelCode::EX_BURIED, datMatrixHolder::EX_BURIED},
    {modelCode::EX_EXPOSED, datMatrixHolder::EX_EXPOSED},
    {modelCode::EHO_EXTENDED, datMatrixHolder::EHO_EXTENDED},
    {modelCode::EHO_HELIX, datMatrixHolder::EHO_HELIX},
    {modelCode::EHO_OTHER, datMatrixHolder::EHO_OTHER},
    {modelCode::EX_EHO_BUR_EXT, datMatrixHolder::EX_EHO_BUR_EXT},
    {modelCode::EX_EHO_BUR_HEL, datMatrixHolder::EX_EHO_BUR_HEL},
    {modelCode::EX_EHO_BUR_OTH, datMatrixHolder::EX_EHO_BUR_OTH},
    {modelCode::EX_EHO_EXP_EXT, datMatrixHolder::EX_EHO_EXP_EXT},
    {modelCode::EX_EHO_EXP_HEL, datMatrixHolder::EX_EHO_EXP_HEL},
    {modelCode::EX_EHO_EXP_OTH, datMatrixHolder::EX_EHO_EXP_OTH}
};


#endif