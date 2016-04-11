//
// Created by Vadim N. on 11/04/2016.
//

#ifndef YMIR_CLONOTYPE_AA_H
#define YMIR_CLONOTYPE_AA_H


#include "clonotype_base.h"
#include "vdj_alignment_aa.h"


namespace ymir {

    class ClonotypeAA : public VDJAlignmentAA, public ClonotypeBase {
    public:


        typedef VDJAlignmentAA vdj_alignment_t;


        ClonotypeAA(const sequence_t &sequence,
                    Recombination recomb,
                    const segments_storage_t &segments,
                    const CodonAlignmentVector &alignments,
                    const n_D_alignments_storage_t &n_D_alignments)
                : ClonotypeBase(sequence,
                                AMINOACID,
                                recomb,
                                segments[0] && segments[1] && ((segments[2] && recomb == VDJ_RECOMB) || recomb == VJ_RECOMB)),
                  VDJAlignmentAA(segments,
                                 alignments,
                                 n_D_alignments)
        {
        }


        ClonotypeAA(const sequence_t &sequence,
                    Recombination recomb,
                    const VDJAlignmentAA &alignment)
                : ClonotypeBase(sequence,
                                AMINOACID,
                                recomb,
                                alignment.nVar() && alignment.nJoi() && ((alignment.nDiv() && recomb == VDJ_RECOMB) || recomb == VJ_RECOMB)),
                  VDJAlignmentAA(alignment)
        {
        }


        virtual ~ClonotypeAA()
        {
        }


    protected:

        ClonotypeAA()
        {
        }

    };
    

    typedef unique_ptr<ClonotypeAA> ClonotypeAAPtr;

}


#endif //YMIR_CLONOTYPE_AA_H
