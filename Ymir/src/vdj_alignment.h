/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vdn at mailbox dot com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _VDJ_ALIGNMENT_H_
#define _VDJ_ALIGNMENT_H_


#include <array>

#include "nogap_alignment_vector.h"


namespace ymir {


    struct VDJAlignment {


        typedef std::array<seg_index_t, 3> segments_storage_t;


        typedef std::vector<seg_index_t> n_D_alignments_storage_t;


        VDJAlignment(const segments_storage_t &segments, 
                     const NoGapAlignmentVector &alignments, 
                     const n_D_alignments_storage_t &n_D_alignments) 
            : _segments(segments), 
              _alignments(alignments),
              _n_D_alignments(n_D_alignments)
        {
        }


        /**
         * \brief Move constructor for _segments, _alignments and _n_D_alignments.
         */
//        VDJAlignment(segments_storage_t &&segments,
//                     NoGapAlignmentVector &&alignments,
//                     n_D_alignments_storage_t &&n_D_alignments)
//            : _segments(std::move(segments)),
//              _alignments(std::move(alignments)),
//              _n_D_alignments(std::move(n_D_alignments))
//        {
//        }


        // VDJAlignment(const VDJAlignment &other) {

        // }


        // VDJAlignment& operator=(const VDJAlignment &other) {

        // }


        virtual ~VDJAlignment()
        {
        }


        bool operator==(const VDJAlignment &other) {
            return _segments == other._segments
                   && _alignments == other._alignments
                   && _n_D_alignments == other._n_D_alignments;
        }

        bool operator!=(const VDJAlignment &other) {
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

        bool isVarMismatch(seg_index_t vgene, seq_len_t pos) const {
            return _alignments.isMismatch(vgene, pos);
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

        bool isJoiMismatch(seg_index_t jgene, seq_len_t pos) const {
            return _alignments.isMismatch(_segments[0] + jgene, pos);
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

        bool isDivMismatch(seg_index_t dgene, seg_index_t align_i, seq_len_t pos) const {
            return _alignments.isMismatch(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i, pos);
        }

        error_num_t numDivMismatches(seg_index_t dgene, seg_index_t align_i, seq_len_t start, seq_len_t end) const {
            return _alignments.numMismatches(_segments[0] + _segments[1] + _n_D_alignments[dgene] + align_i, start, end);
        }
        ///@}


        /**
         * \brief Get the number of D gene alignments for the given D gene.
         */
        seq_len_t numDivAlignments(seg_index_t index) const {
            return _n_D_alignments[index + 1] - _n_D_alignments[index];
        }


    protected:

        segments_storage_t _segments;  /// Two concatenated vectors: vector of length 3 w/ numbers of aligned segments (V-J-D) and
                                              /// vector of indices of segments, aligned on this clone: V1--V2--V3--J1--J2--D1--D2--...

        NoGapAlignmentVector _alignments;  /// Vector of alignments for segments.
        
        n_D_alignments_storage_t _n_D_alignments;  /// Accumulated number of alignments for each D gene segment;
                                                   /// vector's length == _segments[2] + 1. Number of alignments for i-th
                                                   /// gene is equal to v[i+1] - v[i].


        VDJAlignment() {}

    };

}

#endif