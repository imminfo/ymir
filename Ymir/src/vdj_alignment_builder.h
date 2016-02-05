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


    struct VDJAlignmentBuilder : protected VDJAlignment {

        /**
         * \brief Move constructor for _segments, _alignments and _n_D_alignments.
         */
        VDJAlignmentBuilder()
        {
            _segments.fill(0);
        }


        virtual ~VDJAlignmentBuilder()
        {
        }


        VDJAlignment build() {
            segments_storage_t segments;
            segments[0] = _segments[0];
            segments[1] = _segments[1];
            segments[2] = _segments[2];

            n_D_alignments_storage_t nDs(_n_Dalign.size());
            std::copy(_n_Dalign.begin(), _n_Dalign.end(), nDs.begin());

            NoGapAlignmentVector avec;
            avec.extend(_Valign);
            avec.extend(_Jalign);
            avec.extend(_Dalign);

            std::cout << (int) avec.size() << std::endl;

            _Valign.clear();
            _Jalign.clear();
            _Dalign.clear();

            // return std::move(VDJAlignment(std::move(segments), std::move(avec), std::move(nDs)));
            return VDJAlignment(std::move(segments), std::move(avec), std::move(nDs));
        }


        /**
         * \brief Add singular alignments to the builder.
         */
        ///@{
        VDJAlignmentBuilder& addVarAlignment(seg_index_t vseg, seq_len_t vstart, seq_len_t seqstart, seq_len_t alignment_len) {
            _segments[0] += 1;
            _Valign.addAlignment(vseg, vstart, seqstart, alignment_len);
            return *this;
        }


        VDJAlignmentBuilder& addJoiAlignment(seg_index_t jseg, seq_len_t jstart, seq_len_t seqstart, seq_len_t alignment_len) {
            _segments[1] += 1;
            _Jalign.addAlignment(jseg, jstart, seqstart, alignment_len);
            return *this;
        }

        VDJAlignmentBuilder& addDivAlignment(seg_index_t dseg, seq_len_t dstart, seq_len_t seqstart, seq_len_t alignment_len) {
            // check for recombination type?

            if (_Dseg.size() == 0 || dseg != _Dseg[_Dseg.size() - 1]) {
                _n_Dalign.push_back(0);
                _segments[2] += 1;
            }

            _n_Dalign[_n_Dalign.size() - 1] += 1;

            _Dalign.addAlignment(dseg, dstart, seqstart, alignment_len);
            return *this;
        }
        ///@}


        /**
         * \brief Add vectors of alignments to the builder.
         */
        ///@{
        VDJAlignmentBuilder& addVarAlignment(const NoGapAlignmentVector &vec) {
            _Valign.extend(vec);
            return *this;
        }

        VDJAlignmentBuilder& addJoiAlignment(const NoGapAlignmentVector &vec) {
            _Jalign.extend(vec);
            return *this;
        }

        VDJAlignmentBuilder& addDivAlignment(const NoGapAlignmentVector &vec) {
            _Dalign.extend(vec);
            return *this;
        }
        ///@}



    protected:

        NoGapAlignmentVector _Valign, _Jalign, _Dalign;
        std::vector<seg_index_t> _n_Dalign, _Dseg;

    };

}

#endif