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


#include "gapped_alignment_vector.h"
#include "vdj_alignment.h"


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
              _len(events.size())
        {
            _events.swap(events);
        }


        seq_len_t pattern_start() const { return _pattern_start; }


        seq_len_t text_start() const { return _text_start; }


        seq_len_t size() const { return _len; }


    protected:

        seq_len_t _pattern_start, _text_start, _len;
        events_storage_t _events;


        AlignmentBase() {}


        AlignmentBase(seq_len_t p_start, seq_len_t t_start, seq_len_t len) 
            : _pattern_start(p_start), 
              _text_start(t_start), 
              _len(len)
        {
        }

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
            : AlignmentBase(p_start, t_start, len)
        {
        }


        bool isMismatch(seq_len_t i) const { return _events[i]; }


    protected:

        NoGapAlignment() {}

    };


    /**
     *
     */
    struct GappedAlignment : public AlignmentBase {


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


    protected:

        GappedAlignment() {}

    };

}

#endif