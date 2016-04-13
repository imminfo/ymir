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

#ifndef _VDJ_ALIGNMENT_BUILDER_H_
#define _VDJ_ALIGNMENT_BUILDER_H_


#include <iostream>

#include "vdj_alignment.h"


namespace ymir {


    template <typename VDJAlignmentType>
    struct VDJAlignmentBuilder : protected VDJAlignmentBase<typename VDJAlignmentType::alignment_vector_t> {
    public:

        typedef VDJAlignmentType vdj_alignment_t;


        typedef typename vdj_alignment_t::alignment_vector_t alignment_vector_t;


        /**
         * \brief Move constructor for _segments, _alignments and _n_D_alignments.
         */
        VDJAlignmentBuilder()
        {
            vdj_alignment_t::_segments.fill(0);
            _n_Dalign.push_back(0);
        }


        virtual ~VDJAlignmentBuilder()
        {
        }


        VDJAlignmentType buildAlignment() {
            typename VDJAlignmentType::segments_storage_t segments;
            segments[0] = VDJAlignmentType::_segments[0];
            segments[1] = VDJAlignmentType::_segments[1];
            segments[2] = VDJAlignmentType::_segments[2];

            typename VDJAlignmentType::n_D_alignments_storage_t nDs = _n_Dalign;
            for (size_t i = 1; i < nDs.size(); ++i) {
                nDs[i] += nDs[i-1];
            }
            // std::copy(_n_Dalign.begin(), _n_Dalign.end(), nDs.begin());

            alignment_vector_t avec;
            avec.extend(_Valign);
            avec.extend(_Jalign);
            avec.extend(_Dalign);

            _Valign.clear();
            _Jalign.clear();
            _Dalign.clear();

            VDJAlignmentType::_segments.fill(0);
            _n_Dalign.resize(1);

            // return std::move(VDJAlignment(std::move(segments), std::move(avec), std::move(nDs)));
            // return VDJAlignment(std::move(segments), std::move(avec), std::move(nDs));
            return VDJAlignmentType(segments, avec, nDs);
        }


        /**
         * \brief Add singular alignments to the builder.
         */
        ///@{
        VDJAlignmentBuilder& addVarAlignment(seg_index_t vseg, seq_len_t vstart, seq_len_t seqstart, seq_len_t alignment_len) {
            ++vdj_alignment_t::_segments[0];
            _Valign.addAlignment(vseg, vstart, seqstart, alignment_len);
            return *this;
        }


        VDJAlignmentBuilder& addJoiAlignment(seg_index_t jseg, seq_len_t jstart, seq_len_t seqstart, seq_len_t alignment_len) {
            ++vdj_alignment_t::_segments[1];
            _Jalign.addAlignment(jseg, jstart, seqstart, alignment_len);
            return *this;
        }

        VDJAlignmentBuilder& addDivAlignment(seg_index_t dseg, seq_len_t dstart, seq_len_t seqstart, seq_len_t alignment_len) {
            if (_Dseg.size() == 0 || dseg != _Dseg[_Dseg.size() - 1]) {
                _n_Dalign.push_back(0);
                _Dseg.push_back(dseg);
                ++vdj_alignment_t::_segments[2];
            }

            ++_n_Dalign[_n_Dalign.size() - 1];

            _Dalign.addAlignment(dseg, dstart, seqstart, alignment_len);
            return *this;
        }
        ///@}


        /**
         * \brief Add vectors of alignments to the builder.
         *
         * Add alignments to the builder. IMPORTANT NOTE: for V / J alignments each input vectors
         * counts as a set of various V or J gene segments. For D alignments
         * input vec counts as a set of alignments for a specific D gene segment (i.e., only one!).
         */
        ///@{
        VDJAlignmentBuilder& addVarAlignment(const AlignmentVectorBase &vec) {
            _Valign.extend(vec);
            vdj_alignment_t::_segments[0] += vec.size();
            return *this;
        }

        VDJAlignmentBuilder& addJoiAlignment(const AlignmentVectorBase &vec) {
            _Jalign.extend(vec);
            vdj_alignment_t::_segments[1] += vec.size();
            return *this;
        }

        VDJAlignmentBuilder& addDivAlignment(const AlignmentVectorBase &vec) {
            _n_Dalign.push_back(vec.size());
            ++vdj_alignment_t::_segments[2];
            _Dalign.extend(vec);
            return *this;
        }
        ///@}


//        template <GeneSegments GENE>
//        inline VDJAlignmentBuilder& addAlignment(seg_index_t seg_index, seq_len_t genestart, seq_len_t seqstart, seq_len_t alignment_len);

        inline VDJAlignmentBuilder& addAlignment(GeneSegments gene, seg_index_t seg_index, seq_len_t genestart, seq_len_t seqstart, seq_len_t alignment_len) {
            switch (gene) {
                case VARIABLE:
                    return this->addVarAlignment(seg_index, genestart, seqstart, alignment_len);
                case JOINING:
                    return this->addJoiAlignment(seg_index, genestart, seqstart, alignment_len);
                case DIVERSITY:
                    return this->addDivAlignment(seg_index, genestart, seqstart, alignment_len);
                default:
                    return *this;
            }
        }

    protected:

        alignment_vector_t _Valign, _Jalign, _Dalign;
        std::vector<seg_index_t> _n_Dalign, _Dseg;

    };


    typedef VDJAlignmentBuilder<VDJAlignmentNuc> VDJAlignmentNucBuilder;


    typedef VDJAlignmentBuilder<VDJAlignmentAA> VDJAlignmentAABuilder;


//    template <>
//    VDJAlignmentBuilder& VDJAlignmentBuilder::addAlignment<VARIABLE>(seg_index_t seg_index, seq_len_t genestart, seq_len_t seqstart, seq_len_t alignment_len) {
//        return this->addVarAlignment(seg_index, genestart, seqstart, alignment_len);
//    }
//
//
//    template <>
//    VDJAlignmentBuilder& VDJAlignmentBuilder::addAlignment<DIVERSITY>(seg_index_t seg_index, seq_len_t genestart, seq_len_t seqstart, seq_len_t alignment_len) {
//        return this->addDivAlignment(seg_index, genestart, seqstart, alignment_len);
//    }
//
//
//    template <>
//    VDJAlignmentBuilder& VDJAlignmentBuilder::addAlignment<JOINING>(seg_index_t seg_index, seq_len_t genestart, seq_len_t seqstart, seq_len_t alignment_len) {
//        return this->addJoiAlignment(seg_index, genestart, seqstart, alignment_len);
//    }

}

#endif