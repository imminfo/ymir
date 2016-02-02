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


#include "vdj_alignment.h"


namespace ymir {


    struct VDJAlignmentBuilder : protected VDJAlignment {

        /**
         * \brief Move constructor for _segments, _alignments and _n_D_alignments.
         */
        VDJAlignmentBuilder() 
        {
        }


        virtual ~VDJAlignmentBuilder()
        {
        }


        VDJAlignment build() {}


        /**
         * \brief Interface to add alignments to the builder.
         */
        ///@{
        VDJAlignmentBuilder& addVarAlignment(seg_index_t vseg, seq_len_t vstart, seq_len_t seqstart, seq_len_t alignment_len) {
            _Vseg.push_back(vseg);
            _Valign.push_back(vstart);
            _Valign.push_back(seqstart);
            _Valign.push_back(alignment_len);
            return *this;
        }


        VDJAlignmentBuilder& addVarAlignment(const NoGapAlignmentVector &vec) {
            _Valign.extend(vec);
            return *this;
        }

        VDJAlignmentBuilder& addJoiAlignment(seg_index_t jseg, seq_len_t jstart, seq_len_t seqstart, seq_len_t alignment_len) {
            _Jseg.push_back(jseg);
            _Jalign.push_back(jstart);
            _Jalign.push_back(seqstart);
            _Jalign.push_back(alignment_len);
            return *this;
        }

        VDJAlignmentBuilder& addDivAlignment(seg_index_t dseg, seq_len_t dstart, seq_len_t seqstart, seq_len_t alignment_len) {
            // check for recombination type?

            if (_Dseg.size() == 0 || dseg != _Dseg[_Dseg.size() - 1]) {
                _n_Dalign.push_back(0);
            }

            _n_Dalign[_n_Dalign.size() - 1] += 1;

            _Dseg.push_back(dseg);
            _Dalign.push_back(dstart);
            _Dalign.push_back(seqstart);
            _Dalign.push_back(alignment_len);
            return *this;
        }
        ///@}



    protected:


    };

}

#endif