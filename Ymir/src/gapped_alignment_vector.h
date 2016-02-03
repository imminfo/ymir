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

#ifndef _GAPPED_ALIGNMENT_VECTOR_H_
#define _GAPPED_ALIGNMENT_VECTOR_H_


#include "alignment_vector_base.h"


namespace ymir {


    /**
     * \brief A set of functions to modify the vector with various alignment events:
     * matches, mismatches, insertions or deletions.
     */
    ///@{
    inline void add_match(AlignmentVectorBase::events_storage_t *vec)    { vec->push_back(false); vec->push_back(false); }

    inline void add_mismatch(AlignmentVectorBase::events_storage_t *vec) { vec->push_back(false); vec->push_back(true); }

    inline void add_ins(AlignmentVectorBase::events_storage_t *vec)      { vec->push_back(true);  vec->push_back(false); }

    inline void add_del(AlignmentVectorBase::events_storage_t *vec)      { vec->push_back(true);  vec->push_back(true); }
    ///@}


    /**
     *
     */
    struct GappedAlignmentVector : public AlignmentVectorBase {


        GappedAlignmentVector() : AlignmentVectorBase() {

        }


        void addAlignment(seq_len_t p_start, seq_len_t t_start, const events_storage_t &vec, seg_index_t id) {
            _data.push_back(p_start);
            _data.push_back(t_start);
            _data.push_back(vec.size());
            _data.push_back(id);
            _starts.push_back(_events.size());
            _events.insert(_events.end(), vec.begin(), vec.end());
        }


        /**
         * \brief Check if a some alignment has a specific events at a specific position.
         */
        ///@{
        bool isMatch(seq_len_t i, seq_len_t j) const { 
#ifndef DNDEBUG
            if (_starts[i] + j*2 + 1 >= _events.size()) {
                throw(std::runtime_error("Alignment vector: index is out of bounds."));
            }
#endif
            return !(_events[_starts[i] + j*2] && _events[_starts[i] + j*2 + 1]);
        }

        bool isMismatch(seq_len_t i, seq_len_t j) const { 
#ifndef DNDEBUG
            if (_starts[i] + j*2 + 1 >= _events.size()) {
                throw(std::runtime_error("Alignment vector: index is out of bounds."));
            }
#endif
            return !_events[_starts[i] + j*2] && _events[_starts[i] + j*2 + 1];
        }

        bool isIns(seq_len_t i, seq_len_t j) const { 
#ifndef DNDEBUG
            if (_starts[i] + j*2 + 1 >= _events.size()) {
                throw(std::runtime_error("Alignment vector: index is out of bounds."));
            }
#endif
            return _events[_starts[i] + j*2] && !_events[_starts[i] + j*2 + 1];
        }

        bool isDel(seq_len_t i, seq_len_t j) const { 
#ifndef DNDEBUG
            if (_starts[i] + j*2 + 1 >= _events.size()) {
                throw(std::runtime_error("Alignment vector: index is out of bounds."));
            }
#endif
            return _events[_starts[i] + j*2] && _events[_starts[i] + j*2 + 1];
        }
        ///@}

    };

}

#endif