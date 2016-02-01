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


    struct NoGapAlignment {

        typedef std::vector<bool> mismatch_storage_t;


        NoGapAlignment(seq_len_t p_start, seq_len_t t_start, mismatch_storage_t &errors) 
            : _pattern_start(p_start), 
              _text_start(t_start), 
              _len(errors.size())
        {
            _errors.swap(errors);
        }


        NoGapAlignment(seq_len_t p_start, seq_len_t t_start, seq_len_t len) 
            : _pattern_start(p_start), 
              _text_start(t_start), 
              _len(len)
        {
        }


        seq_len_t pattern_start() const { return _pattern_start; }


        seq_len_t text_start() const { return _text_start; }


        seq_len_t len() const { return _len; }


        bool isMismatch(seq_len_t i) const { return _errors[i]; }


    private:

        seq_len_t _pattern_start, _text_start, _len;
        mismatch_storage_t _errors;


        NoGapAlignment() {}

    };


    struct NoGapAlignmentVector {

        // reserve

        seq_len_t pattern_start(seq_len_t i) const {}

        seq_len_t text_start(seq_len_t i) const {}

        seq_len_t len(seq_len_t i) const { return _len; }

        bool isMismatch(seq_len_t i, seq_len_t j) const { return _errors[j]; }

        void addAlignment() {}

        void finish();

    private:


    };


    struct GappedAlignment {

        /**
         * \typedef alignment_events_storage_t
         */
        typedef std::vector<bool> alignment_events_storage_t;


        GappedAlignment(seq_len_t p_start, seq_len_t t_start, alignment_events_storage_t &events) 
            : _pattern_start(p_start), 
              _text_start(t_start)
        {
            _events.swap(events);
        }


        seq_len_t pattern_start() const { return _pattern_start; }


        seq_len_t text_start() const { return _text_start; }


        ///@{
        bool isMatch(seq_len_t i) const { return !(_events[i*2] && _events[i*2 + 1]); }

        bool isMismatch(seq_len_t i) const { return !_events[i*2] && _events[i*2 + 1]; }

        bool isIns(seq_len_t i) const { return _events[i*2] && !_events[i*2 + 1]; }

        bool isDel(seq_len_t i) const { return _events[i*2] && _events[i*2 + 1]; }
        ///@}


    private:

        seq_len_t _pattern_start, _text_start;
        Events _events;


        GappedAlignment() {}

    };



    ///@{
    inline void add_match(GappedAlignment::alignment_events_storage_t *vec)    { vec->push_back(false); vec->push_back(false); }

    inline void add_mismatch(GappedAlignment::alignment_events_storage_t *vec) { vec->push_back(false); vec->push_back(true); }

    inline void add_ins(GappedAlignment::alignment_events_storage_t *vec)      { vec->push_back(true);  vec->push_back(false); }

    inline void add_del(GappedAlignment::alignment_events_storage_t *vec)      { vec->push_back(true);  vec->push_back(true); }
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