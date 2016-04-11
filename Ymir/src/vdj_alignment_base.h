//
// Created by Vadim N. on 11/04/2016.
//

#ifndef YMIR_VDJ_ALIGNMENT_BASE_H
#define YMIR_VDJ_ALIGNMENT_BASE_H


#include <array>


namespace ymir {

    template <typename AlignmentVectorType>
    class VDJAlignmentBase {
    public:

        typedef std::array<seg_index_t, 3> segments_storage_t;


        typedef std::vector<seg_index_t> n_D_alignments_storage_t;


        VDJAlignmentBase(const segments_storage_t &segments,
                         const AlignmentVectorType &alignments,
                         const n_D_alignments_storage_t &n_D_alignments)
            : _segments(segments),
              _alignments(alignments),
              _n_D_alignments(n_D_alignments)
        {
        }


        virtual ~VDJAlignmentBase()
        {
        }


        bool operator==(const VDJAlignmentBase &other) {
            return _segments == other._segments
                   && _alignments == other._alignments
                   && _n_D_alignments == other._n_D_alignments;
        }

        bool operator!=(const VDJAlignmentBase &other) {
            return _segments != other._segments
                   || _alignments != other._alignments
                   || _n_D_alignments != other._n_D_alignments;
        }


        /**
         * \brief Get the number of alignments for the specific gene.
         */
        ///@{
        seg_index_t nVar() const { return _segments[0]; }

        seg_index_t nJoi() const { return _segments[1]; }

        seg_index_t nDiv() const { return _segments[2]; }
        ///@}


        /**
         * \brief Get the index of the aligned gene segment for the specific gene.
         */
        ///@{
        seg_index_t getVar(size_t index) const { return _alignments.id(index); }

        seg_index_t getJoi(size_t index) const { return _alignments.id(_segments[0] + index); }

        seg_index_t getDiv(size_t index) const { return _alignments.id(_segments[0] + _segments[1] + _n_D_alignments[index]); }
        ///@}


        /**
         * \brief Get alignments for the specific gene.
         */
        ///@{
        seq_len_t getVarGeneStart(seg_index_t vgene) const {
            return _alignments.pattern_start(vgene);
        }

        seq_len_t getVarGeneEnd(seg_index_t vgene) const {
            return _alignments.pattern_start(vgene) + _alignments.len(vgene) - 1;
        }

        seq_len_t getVarSeqStart(seg_index_t vgene) const {
            return _alignments.text_start(vgene);
        }

        seq_len_t getVarSeqEnd(seg_index_t vgene) const {
            return _alignments.text_start(vgene) + _alignments.len(vgene) - 1;
        }

        seq_len_t getVarLen(seg_index_t vgene) const {
            return _alignments.len(vgene);
        }


        seq_len_t getJoiGeneStart(seg_index_t jgene) const {
            return _alignments.pattern_start(_segments[0] + jgene);
        }

        seq_len_t getJoiGeneEnd(seg_index_t jgene) const {
            return _alignments.pattern_start(_segments[0] + jgene) + _alignments.len(_segments[0] + jgene) - 1;
        }

        seq_len_t getJoiSeqStart(seg_index_t jgene) const {
            return _alignments.text_start(_segments[0] + jgene);
        }

        seq_len_t getJoiSeqEnd(seg_index_t jgene) const {
            return _alignments.text_start(_segments[0] + jgene) + _alignments.len(_segments[0] + jgene) - 1;
        }

        seq_len_t getJoiLen(seg_index_t jgene) const {
            return _alignments.len(_segments[0] + jgene);
        }


        seq_len_t getDivGeneStart(seg_index_t dgene, seg_index_t align_i) const {
            return _alignments.pattern_start(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i);
        }

        seq_len_t getDivGeneEnd(seg_index_t dgene, seg_index_t align_i) const {
            return _alignments.pattern_start(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i)
                   + _alignments.len(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i) - 1;
        }

        seq_len_t getDivSeqStart(seg_index_t dgene, seg_index_t align_i) const {
            return _alignments.text_start(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i);
        }

        seq_len_t getDivSeqEnd(seg_index_t dgene, seg_index_t align_i) const {
            return _alignments.text_start(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i)
                   + _alignments.len(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i) - 1;
        }

        seq_len_t getDivLen(seg_index_t dgene, seg_index_t align_i) const {
            return _alignments.len(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i);
        }
        ///@}


        /**
         * \brief Get the number of D gene alignments for the given D gene.
         */
        seq_len_t numDivAlignments(seg_index_t index) const {
            return _n_D_alignments[index + 1] - _n_D_alignments[index];
        }


    protected:

        AlignmentVectorType _alignments;  /// Vector of alignments for segments.

        n_D_alignments_storage_t _n_D_alignments;  /// Accumulated number of alignments for each D gene segment;
                                                   /// vector's length == _segments[2] + 1. Number of alignments for i-th
                                                   /// gene is equal to v[i+1] - v[i].

        segments_storage_t _segments;  /// Two concatenated vectors: vector of length 3 w/ numbers of aligned segments (V-J-D) and
                                       /// vector of indices of segments, aligned on this clone: V1--V2--V3--J1--J2--D1--D2--...


        VDJAlignmentBase() {}

    };
}


#endif //YMIR_VDJ_ALIGNMENT_BASE_H
