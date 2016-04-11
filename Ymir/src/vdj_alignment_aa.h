//
// Created by Vadim N. on 11/04/2016.
//

#ifndef YMIR_VDJ_ALIGNMENT_AA_H
#define YMIR_VDJ_ALIGNMENT_AA_H


#include "codon_alignment_vector.h"

#include "vdj_alignment_base.h"


namespace ymir {

    class VDJAlignmentAA : public VDJAlignmentBase<CodonAlignmentVector> {

    public:

        typedef VDJAlignmentAA base_t;


        VDJAlignmentAA(const segments_storage_t &segments,
                       const CodonAlignmentVector &alignments,
                       const n_D_alignments_storage_t &n_D_alignments)
            : VDJAlignmentBase<CodonAlignmentVector>(segments, alignments, n_D_alignments)
        {
        }


        virtual ~VDJAlignmentAA()
        {
        }


        //
        // Working with codons is here
        //


    protected:

        VDJAlignmentAA() : VDJAlignmentBase()
        {
        }

    };

}

#endif //YMIR_VDJ_ALIGNMENT_AA_H
