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

        static const SequenceType sequence_type = AMINOACID;


        typedef VDJAlignmentAA vdj_alignment_t;


        ClonotypeAA(const sequence_t &sequence,
                    Recombination recomb,
                    const segments_storage_t &segments,
                    const CodonAlignmentVector &alignments,
                    const n_D_alignments_storage_t &n_D_alignments)
                : ClonotypeBase(sequence,
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
                                recomb,
                                alignment.nVar() && alignment.nJoi() && ((alignment.nDiv() && recomb == VDJ_RECOMB) || recomb == VJ_RECOMB)),
                  VDJAlignmentAA(alignment)
        {
        }


        virtual ~ClonotypeAA()
        {
        }


        /**
         * \brief Check if clonotype's sequence is coding, noncoding or out-of-frame.
         */
        ///@{
        bool isCoding() const { return !has_bad_aa_codons(_sequence); }

        bool isNoncoding() const { return has_bad_aa_codons(_sequence); }

        bool isOutOfFrame() const { return has_oof_aa_codon(_sequence); }
        ///@}


    protected:

        ClonotypeAA()
        {
        }

    };
    

    typedef unique_ptr<ClonotypeAA> ClonotypeAAPtr;

}


#endif //YMIR_CLONOTYPE_AA_H
