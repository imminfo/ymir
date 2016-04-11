//
// Created by Vadim N. on 11/04/2016.
//

#ifndef YMIR_CLONOTYPE_NUC_H
#define YMIR_CLONOTYPE_NUC_H


#include "clonotype_base.h"
#include "vdj_alignment_nuc.h"


namespace ymir {


    class ClonotypeNuc : public VDJAlignmentNuc, ClonotypeBase {
    public:

        ClonotypeNuc(const sequence_t &sequence,
                     Recombination recomb,
                     const segments_storage_t &segments,
                     const NoGapAlignmentVector &alignments,
                     const n_D_alignments_storage_t &n_D_alignments)
            : ClonotypeBase(sequence,
                            NUCLEOTIDE,
                            recomb,
                            segments[0] && segments[1] && ((segments[2] && recomb == VDJ_RECOMB) || recomb == VJ_RECOMB)),
              VDJAlignmentNuc(segments,
                              alignments,
                              n_D_alignments)
        {
        }


        ClonotypeNuc(const sequence_t &sequence,
                     Recombination recomb,
                     const VDJAlignmentNuc &alignment)
            : ClonotypeBase(sequence,
                            NUCLEOTIDE,
                            recomb,
                            alignment.nVar() && alignment.nJoi() && ((alignment.nDiv() && recomb == VDJ_RECOMB) || recomb == VJ_RECOMB)),
              VDJAlignmentNuc(alignment)
        {
        }


        virtual ~ClonotypeNuc()
        {
        }


    protected:

        ClonotypeNuc()
        {
        }

    };


    typedef unique_ptr<ClonotypeNuc> ClonotypeNucPtr;
}


#endif //YMIR_CLONOTYPE_NUC_H
