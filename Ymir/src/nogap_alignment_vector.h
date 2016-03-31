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

#ifndef _NOGAP_ALIGNMENT_VECTOR_H_
#define _NOGAP_ALIGNMENT_VECTOR_H_


#include "alignment_vector_base.h"
#include "multimatrixchain.h"


namespace ymir {

    /**
     *
     */
    struct NoGapAlignmentVector : public AlignmentVectorBase {


        NoGapAlignmentVector() : AlignmentVectorBase()
        {
        }


        /**
         * \brief Add a new alignment to the vector.
         * 
         * \param id
         * \param p_start
         * \param t_start
         * \param size
         * \param vec
         */
        ///@{
        void addAlignment(seg_index_t id, seq_len_t p_start, seq_len_t t_start, seq_len_t size) {
            _data.push_back(p_start);
            _data.push_back(t_start);
            _data.push_back(size);
            _data.push_back(id);
        }

        void addAlignment(seg_index_t id, seq_len_t p_start, seq_len_t t_start, const events_storage_t &vec) {
            _data.push_back(p_start);
            _data.push_back(t_start);
            _data.push_back(vec.size());
            _data.push_back(id);
            _starts.push_back(_events.size());
            _events.insert(_events.end(), vec.begin(), vec.end());
        }
        ///@}


        /**
         *
         */
        bool isMismatch(seq_len_t i, seq_len_t j) const {
#ifndef DNDEBUG
            check_and_throw(_events.size() == 0, "Alignment vector: no errors stored in the vector.");

            if (_starts[i] + j - 1 >= _events.size()) {
                std::cout << (int) _starts.size() << std::endl;
                std::cout << (int) i << std::endl;
                std::cout << (int) _starts[i] << std::endl;
                std::cout << (int) j << std::endl;
                std::cout << (int) _events.size() << std::endl;
                throw(std::runtime_error("Alignment vector: mismatch index is out of bounds."));
            }
#endif
            return _events[_starts[i] + j - 1];
        }


        error_num_t numMismatches(seq_len_t i, seq_len_t start, seq_len_t end) const {
#ifndef DNDEBUG
            check_and_throw(_events.size() == 0, "Alignment vector: no errors stored in the vector.");

            if (_starts[i] + start - 1 >= _events.size()) {
                throw(std::runtime_error("Alignment vector: start index is out of bounds."));
            }
            if (_starts[i] + end - 1 >= _events.size()) {
                throw(std::runtime_error("Alignment vector: end index is out of bounds."));
            }
#endif
            error_num_t res = 0;
            for (size_t j = start; j <= end; ++j) {
                res += static_cast<error_num_t>(_events[_starts[i] + j - 1]);
            }
            return res;
        }

    };

}

#endif