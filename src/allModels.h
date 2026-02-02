
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

#endif