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

#include "tools.h"
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


        // waiting until C++14 for full support?
        // AlignmentVectorBase(AlignmentVectorBase &&other) 
        //     : _data(std::move(other._data)),
        //       _events(std::move(other._events)), 
        //       _starts(std::move(other._starts))
        // {
        // }

        bool operator==(const AlignmentVectorBase &other) {
            return _data == other._data
                   && _events == other._events
                   && _starts == other._starts;
        }

        bool operator!=(const AlignmentVectorBase &other) {
            return _data != other._data
                   || _events != other._events
                   || _starts != other._starts;
        }


        void extend(const AlignmentVectorBase &other) {
            if (other.size()) {
                _data.reserve(_data.size() + other._data.size() + 1);
                _data.insert(_data.end(), other._data.begin(), other._data.end());

                _starts.reserve(_starts.size() + other._starts.size() + 1);
                size_t events_size = _events.size();
                for (size_t i = 0; i < other._starts.size(); ++i) {
                    _starts.push_back(events_size + other._starts[i]);
                }

                _events.reserve(_events.size() + other._events.size() + 1);
                _events.insert(_events.end(), other._events.begin(), other._events.end());
            }
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
            check_and_throw(i*4 >= _data.size(), "Alignment vector: pattern index " + std::to_string(i*4) + " is out of bounds (" + std::to_string(_data.size()) + ")");
            return _data[i*4];
        }


        seq_len_t text_start(seq_len_t i) const {
            check_and_throw(i*4 + 1 >= _data.size(), "Alignment vector: text index " + std::to_string(i*4 + 1) + " is out of bounds (" + std::to_string(_data.size()) + ")");
            return _data[i*4 + 1];
        }


        seq_len_t len(seq_len_t i) const {
            check_and_throw(i*4 + 2 >= _data.size(), "Alignment vector: length index " + std::to_string(i*4 + 2) + " is out of bounds (" + std::to_string(_data.size()) + ")");
            return _data[i*4 + 2];
        }


        seq_len_t id(seq_len_t i) const {
            check_and_throw(i*4 + 3 >= _data.size(), "Alignment vector: ID index " + std::to_string(i*4 + 3) + " is out of bounds (" + std::to_string(_data.size()) + ")");
            return _data[i*4 + 3];
        }


        bool hasEvents() const { return static_cast<bool>(_events.size()); }


    protected:

        std::vector<seq_len_t> _data;  /// Vector of 4-tpuples - pattern start, text start, alignment length and a text ID.
        events_storage_t _events;
        std::vector<size_t> _starts;

    };

}

#endif