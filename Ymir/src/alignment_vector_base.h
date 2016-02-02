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

#ifndef _ALIGNMENT_VECTOR_BASE_H_
#define _ALIGNMENT_VECTOR_BASE_H_


#include <string>
#include <vector>

#include "types.h"


namespace ymir {

    /**
     *
     */
    struct AlignmentVectorBase {


        typedef std::vector<bool> events_storage_t;

            
        static const size_t default_start_reserve_size = 20;

        static const size_t default_events_reserve_size = 600;


        AlignmentVectorBase() {
            _events.reserve(default_events_reserve_size);
            _starts.reserve(default_start_reserve_size);
        }


        // waiting until C++14 ?
        AlignmentVectorBase(AlignmentVectorBase &&other) 
            : _events(std::move(other._events)), 
              _starts(std::move(other._starts))
        {
        }


        void extend(const AlignmentVectorBase &other) {
            _data.insert(_data.begin(), other._data.begin(), other._data.end());
            _starts.insert(_starts.begin(), other._starts.begin(), other._starts.end());
            _events.insert(_events.begin(), other._events.begin(), other._events.end());
        }


        size_t size() const { 
            return _data.size() / 3;
        }


        seq_len_t pattern_start(seq_len_t i) const { 
#ifndef DNDEBUG
            if (i*3 >= _data.size()) {
                throw(std::runtime_error("Alignment vector: index is out of bounds."));
            }
#endif
            return _data[i*3];
        }


        seq_len_t text_start(seq_len_t i) const {
#ifndef DNDEBUG
            if (i*3 + 1 >= _data.size()) {
                throw(std::runtime_error("Alignment vector: index is out of bounds."));
            }
#endif
            return _data[i*3 + 1];
        }


        seq_len_t len(seq_len_t i) const { 
#ifndef DNDEBUG
            if (i*3 + 2 >= _data.size()) {
                throw(std::runtime_error("Alignment vector: index is out of bounds."));
            }
#endif
            return _data[i*3 + 2];
        }


        void finish() {
            _data.reserve(_data.size() + 1);
            _events.reserve(_events.size() + 1);
        }


    protected:

        std::vector<seq_len_t> _data;
        events_storage_t _events;;
        std::vector<size_t> _starts;

    };

}

#endif