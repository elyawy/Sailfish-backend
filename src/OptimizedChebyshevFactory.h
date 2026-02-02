// OptimizedChebyshevFactory.h
// Helper to create chebyshevAccelerator with precomputed data when available

#ifndef OPTIMIZED_CHEBYSHEV_FACTORY_H
#define OPTIMIZED_CHEBYSHEV_FACTORY_H

#include "../libs/Phylolib/includes/chebyshevAccelerator.h"
#include "../libs/Phylolib/includes/replacementModel.h"
#include "../libs/Phylolib/includes/PrecomputedChebyshevData.h"
#include "allModels.h"
#include <memory>

class OptimizedChebyshevFactory {
public:
    // Create a chebyshevAccelerator, using precomputed data if available
    static std::unique_ptr<chebyshevAccelerator> create(
        replacementModel* model, 
        modelCode code
    ) {
        // Check if we have precomputed data for this model
        switch (code) {
            case CPREV45:
                return createFromPrecomputed(model, 
                    PrecomputedCheby::cpREV45_coff,
                    PrecomputedCheby::cpREV45_derv_coff,
                    PrecomputedCheby::cpREV45_sec_derv_coff);
            
            case DAYHOFF:
                return createFromPrecomputed(model,
                    PrecomputedCheby::dayhoff_coff,
                    PrecomputedCheby::dayhoff_derv_coff,
                    PrecomputedCheby::dayhoff_sec_derv_coff);
            
            case JONES: // JTT
                return createFromPrecomputed(model,
                    PrecomputedCheby::jones_coff,
                    PrecomputedCheby::jones_derv_coff,
                    PrecomputedCheby::jones_sec_derv_coff);
            
            case MTREV24:
                return createFromPrecomputed(model,
                    PrecomputedCheby::mtREV24_coff,
                    PrecomputedCheby::mtREV24_derv_coff,
                    PrecomputedCheby::mtREV24_sec_derv_coff);
            
            case WAG:
                return createFromPrecomputed(model,
                    PrecomputedCheby::wag_coff,
                    PrecomputedCheby::wag_derv_coff,
                    PrecomputedCheby::wag_sec_derv_coff);
            
            case HIVB:
                return createFromPrecomputed(model,
                    PrecomputedCheby::HIVb_coff,
                    PrecomputedCheby::HIVb_derv_coff,
                    PrecomputedCheby::HIVb_sec_derv_coff);
            
            case HIVW:
                return createFromPrecomputed(model,
                    PrecomputedCheby::HIVw_coff,
                    PrecomputedCheby::HIVw_derv_coff,
                    PrecomputedCheby::HIVw_sec_derv_coff);
            
            case LG:
                return createFromPrecomputed(model, 
                    PrecomputedCheby::lg_coff,
                    PrecomputedCheby::lg_derv_coff,
                    PrecomputedCheby::lg_sec_derv_coff);
            
            // NOTE: EMPIRICODON skipped (codon model, uses normal constructor)
            
            case EX_BURIED:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EX_BURIED_coff,
                    PrecomputedCheby::EX_BURIED_derv_coff,
                    PrecomputedCheby::EX_BURIED_sec_derv_coff);
            
            case EX_EXPOSED:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EX_EXPOSED_coff,
                    PrecomputedCheby::EX_EXPOSED_derv_coff,
                    PrecomputedCheby::EX_EXPOSED_sec_derv_coff);
            
            case EHO_EXTENDED:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EHO_EXTENDED_coff,
                    PrecomputedCheby::EHO_EXTENDED_derv_coff,
                    PrecomputedCheby::EHO_EXTENDED_sec_derv_coff);
            
            case EHO_HELIX:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EHO_HELIX_coff,
                    PrecomputedCheby::EHO_HELIX_derv_coff,
                    PrecomputedCheby::EHO_HELIX_sec_derv_coff);
            
            case EHO_OTHER:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EHO_OTHER_coff,
                    PrecomputedCheby::EHO_OTHER_derv_coff,
                    PrecomputedCheby::EHO_OTHER_sec_derv_coff);
            
            case EX_EHO_BUR_EXT:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EX_EHO_BUR_EXT_coff,
                    PrecomputedCheby::EX_EHO_BUR_EXT_derv_coff,
                    PrecomputedCheby::EX_EHO_BUR_EXT_sec_derv_coff);
            
            case EX_EHO_BUR_HEL:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EX_EHO_BUR_HEL_coff,
                    PrecomputedCheby::EX_EHO_BUR_HEL_derv_coff,
                    PrecomputedCheby::EX_EHO_BUR_HEL_sec_derv_coff);
            
            case EX_EHO_BUR_OTH:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EX_EHO_BUR_OTH_coff,
                    PrecomputedCheby::EX_EHO_BUR_OTH_derv_coff,
                    PrecomputedCheby::EX_EHO_BUR_OTH_sec_derv_coff);
            
            case EX_EHO_EXP_EXT:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EX_EHO_EXP_EXT_coff,
                    PrecomputedCheby::EX_EHO_EXP_EXT_derv_coff,
                    PrecomputedCheby::EX_EHO_EXP_EXT_sec_derv_coff);
            
            case EX_EHO_EXP_HEL:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EX_EHO_EXP_HEL_coff,
                    PrecomputedCheby::EX_EHO_EXP_HEL_derv_coff,
                    PrecomputedCheby::EX_EHO_EXP_HEL_sec_derv_coff);
            
            case EX_EHO_EXP_OTH:
                return createFromPrecomputed(model,
                    PrecomputedCheby::EX_EHO_EXP_OTH_coff,
                    PrecomputedCheby::EX_EHO_EXP_OTH_derv_coff,
                    PrecomputedCheby::EX_EHO_EXP_OTH_sec_derv_coff);
                
            default:
                // Custom or other models - use normal constructor (30-60ms computation)
                return std::make_unique<chebyshevAccelerator>(model);
        }
    }

private:
    // Create using precomputed coefficients (instant, no computation)
    static std::unique_ptr<chebyshevAccelerator> createFromPrecomputed(
        replacementModel* model,
        const VVVdouble& coff,
        const VVVdouble& derv_coff,
        const VVVdouble& sec_derv_coff
    ) {
        return std::make_unique<chebyshevAccelerator>(
            model,
            coff,
            derv_coff,
            sec_derv_coff
            // Default parameters: alphanetSize=20, totalNumOfCoef=60, 
            // usingNumberOfCoef=13, rightRange=0, leftRange=2
        );
    }
};

#endif // OPTIMIZED_CHEBYSHEV_FACTORY_H