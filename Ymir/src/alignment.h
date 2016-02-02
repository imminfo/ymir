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

#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_


#include <string>
#include <vector>

#include "types.h"


namespace ymir {

    #define DEFAULT_LOCAL_ALIGNMENT_RESERVE 160


    /**
     * \typedef alignment_score_t
     *
     * \brief Type for scores of gene segments alignments to input sequences.
     */
    typedef int16_t alignment_score_t;


    /**
     *
     */
    struct AlignmentBase {

        /**
         * \typedef events_storage_t
         */
        typedef std::vector<bool> events_storage_t;


        AlignmentBase(seq_len_t p_start, seq_len_t t_start, events_storage_t &events) 
            : _pattern_start(p_start), 
              _text_start(t_start), 
              _len(errors.size())
        {
            _errors.swap(errors);
        }


        seq_len_t pattern_start() const { return _pattern_start; }


        seq_len_t text_start() const { return _text_start; }


        seq_len_t len() const { return _len; }


    private:

        seq_len_t _pattern_start, _text_start, _len;
        events_storage_t _events;


        AlignmentBase() {}

    };


    /**
     *
     */
    struct NoGapAlignment : public AlignmentBase {

        /**
         *
         */
        NoGapAlignment(seq_len_t p_start, seq_len_t t_start, events_storage_t &events) 
            : AlignmentBase(p_start, t_start, events)
        {
        }


        NoGapAlignment(seq_len_t p_start, seq_len_t t_start, seq_len_t len) 
            : _pattern_start(p_start), 
              _text_start(t_start), 
              _len(len)
        {
        }


        bool isMismatch(seq_len_t i) const { return _errors[i]; }


    private:

        NoGapAlignment() {}

    };


    /**
     *
     */
    struct GappedAlignment {


        GappedAlignment(seq_len_t p_start, seq_len_t t_start, events_storage_t &events) 
            : AlignmentBase(p_start, t_start, events)
        {
        }


        ///@{
        bool isMatch(seq_len_t i) const { return !(_events[i*2] && _events[i*2 + 1]); }

        bool isMismatch(seq_len_t i) const { return !_events[i*2] && _events[i*2 + 1]; }

        bool isIns(seq_len_t i) const { return _events[i*2] && !_events[i*2 + 1]; }

        bool isDel(seq_len_t i) const { return _events[i*2] && _events[i*2 + 1]; }
        ///@}


    private:

        GappedAlignment() {}

    };


    /**
     *
     */
    struct AlignmentVectorBase {


        typedef std::vector<bool> events_storage_t;

            
        static const size_t default_reserve_size = 100;


        AlignmentVectorBase() {

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


        bool isMismatch(seq_len_t i, seq_len_t j) const { 
#ifndef DNDEBUG
            if (_starts[i] + j >= _data.size()) {
                throw(std::runtime_error("Alignment vector: index is out of bounds."));
            }
#endif
            return _errors[_starts[i] + j];
        }


        void addAlignment(const mismatch_storage_t &vec) {
            _starts.push_back(_errors.size());
            _errors.insert(_errors.end(), vec.begin(), vec.end());
        }


        void finish() {
            _data.reserve(_data.size() + 1);
            _errors.reserve(_errors.size() + 1);
        }


    private:

        std::vector<seq_len_t> _data;
        events_storage_t _events;;
        std::vector<seq_len_t> _starts;

    };


    ///@{
    inline void add_match(GappedAlignment::events_storage_t *vec)    { vec->push_back(false); vec->push_back(false); }

    inline void add_mismatch(GappedAlignment::events_storage_t *vec) { vec->push_back(false); vec->push_back(true); }

    inline void add_ins(GappedAlignment::events_storage_t *vec)      { vec->push_back(true);  vec->push_back(false); }

    inline void add_del(GappedAlignment::events_storage_t *vec)      { vec->push_back(true);  vec->push_back(true); }
    ///@}


    /**
     * \struct Alignment
     *
     */
    struct SegmentAlignment {

        seg_index_t segment;
        seq_len_t start, end;
        alignment_score_t score;
        AlignmentEventVector events;


        SegmentAlignment()
                : segment(0), start(0), end(0), score(0)
        { }


        SegmentAlignment(seg_index_t segment_, seq_len_t start_, seq_len_t end_, alignment_score_t score_, const AlignmentEventVector &events_)
                : segment(segment_), start(start_), end(end_), score(score_), events(events_)
        { }

    };


    typedef std::vector<SegmentAlignment> SegmentAlignmentVector;

}

#endif