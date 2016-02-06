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
            _starts.push_back(0);
        }


        // waiting until C++14 for full support?
        // AlignmentVectorBase(AlignmentVectorBase &&other) 
        //     : _data(std::move(other._data)),
        //       _events(std::move(other._events)), 
        //       _starts(std::move(other._starts))
        // {
        // }


        void extend(const AlignmentVectorBase &other) {
            _data.reserve(_data.size() + other._data.size() + 1);
            _data.insert(_data.end(), other._data.begin(), other._data.end());

            _starts.reserve(_starts.size() + other._starts.size() + 1);
            _starts.insert(_starts.end(), other._starts.begin(), other._starts.end());

            _events.reserve(_events.size() + other._events.size() + 1);
            for (size_t i = 0; i < other._events.size(); ++i) {
                _events.push_back(other._events[i]);
            }
            // _events.insert(_events.end(), other._events.begin(), other._events.end());
        }


        size_t size() const { 
            return _data.size() / 4;
        }


        void finish() {
            _data.reserve(_data.size() + 1);
            _events.reserve(_events.size() + 1);
        }


        void clear() {
            _data.clear();
            _starts.clear();
            _events.clear();
        }


        seq_len_t pattern_start(seq_len_t i) const { 
#ifndef DNDEBUG
            if (i*4 >= _data.size()) {
                throw(std::runtime_error("Alignment vector: pattern index is out of bounds."));
            }
#endif
            return _data[i*4];
        }


        seq_len_t text_start(seq_len_t i) const {
#ifndef DNDEBUG
            if (i*4 + 1 >= _data.size()) {
                throw(std::runtime_error("Alignment vector: test index is out of bounds."));
            }
#endif
            return _data[i*4 + 1];
        }


        seq_len_t len(seq_len_t i) const { 
#ifndef DNDEBUG
            if (i*4 + 2 >= _data.size()) {
                throw(std::runtime_error("Alignment vector: length index is out of bounds."));
            }
#endif
            return _data[i*4 + 2];
        }


        seq_len_t id(seq_len_t i) const { 
#ifndef DNDEBUG
            if (i*4 + 3 >= _data.size()) {
                throw(std::runtime_error("Alignment vector: ID index is out of bounds."));
            }
#endif
            return _data[i*4 + 3];
        }


    protected:

        std::vector<seq_len_t> _data;  /// Vector of 4-tpuples - pattern start, text start, alignment length and a text ID.
        events_storage_t _events;
        std::vector<size_t> _starts;

    };

}

#endif