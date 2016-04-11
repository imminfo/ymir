//
// Created by Vadim N. on 11/04/2016.
//

#ifndef YMIR_VDJ_ALIGNMENT_NUC_H
#define YMIR_VDJ_ALIGNMENT_NUC_H


#include "nogap_alignment_vector.h"

#include "vdj_alignment_base.h"


namespace ymir {


    class VDJAlignmentNuc : public VDJAlignmentBase<NoGapAlignmentVector> {

    public:

        VDJAlignmentNuc(const segments_storage_t &segments,
                        const NoGapAlignmentVector &alignments,
                        const n_D_alignments_storage_t &n_D_alignments)
            : VDJAlignmentBase<NoGapAlignmentVector>(segments, alignments, n_D_alignments)
        {
        }


        virtual ~VDJAlignmentNuc()
        {
        }


        ///@{
        bool isVarMismatch(seg_index_t vgene, seq_len_t pos) const {
            return _alignments.isMismatch(vgene, pos);
        }

        bool isJoiMismatch(seg_index_t jgene, seq_len_t pos) const {
            return _alignments.isMismatch(_segments[0] + jgene, pos);
        }

        bool isDivMismatch(seg_index_t dgene, seg_index_t align_i, seq_len_t pos) const {
            return _alignments.isMismatch(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i, pos);
        }

        error_num_t numDivMismatches(seg_index_t dgene, seg_index_t align_i, seq_len_t start, seq_len_t end) const {
            return _alignments.numMismatches(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i, start, end);
        }
        ///@}


    protected:

        VDJAlignmentNuc() : VDJAlignmentBase()
        {
        }

    };
}

#endif //YMIR_VDJ_ALIGNMENT_NUC_H
