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


        VDJAlignmentAA(const segments_storage_t &segments,
                       const CodonAlignmentVector &alignments,
                       const n_D_alignments_storage_t &n_D_alignments)
            : VDJAlignmentBase<CodonAlignmentVector>(segments, alignments, n_D_alignments)
        {
        }


        virtual ~VDJAlignmentAA()
        {
        }


        /**
         *
         */
        ///@{
        codon_hash getVarCodon(seg_index_t vgene, seq_len_t pos) const {
            return _alignments.getCodon(vgene, pos);
        }

        codon_hash getDivCodon(seg_index_t dgene, seg_index_t align_i, seq_len_t pos) const {
            return _alignments.getCodon(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i, pos);
        }

        codon_hash getJoiCodon(seg_index_t jgene, seq_len_t pos) const {
            return _alignments.getCodon(_segments[0] + jgene, pos);
        }
        ///@}


    protected:

        VDJAlignmentAA() : VDJAlignmentBase()
        {
        }

    };

}

#endif //YMIR_VDJ_ALIGNMENT_AA_H
