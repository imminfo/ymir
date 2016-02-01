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

        typedef std::vector<bool> Mismatches;


        NoGapAlignment(seq_len_t p_start, seq_len_t t_start, const Mismatches &errors) 
            : _pattern_start(p_start), 
              _text_start(t_start), 
              _len(errors.size()),
              _errors(errors)
        {
        }


        NoGapAlignment(seq_len_t p_start, seq_len_t t_start, seq_len_t len) 
            : _pattern_start(p_start), 
              _text_start(t_start), 
              _len(len)
        {
        }


        seq_len_t pattern_start() const { return _pattern_start; }


        seq_len_t _text_start() const { return _text_start; }


        seq_len_t len() const { return _len; }


        bool isMismatch(seq_len_t i) const { return _errors[i]; }


    private:

        seq_len_t _pattern_start, _text_start, _len;
        Mismatches _errors;


        NoGapAlignment() {}

    };


    struct AlignmentEvents {

    };


    struct GappedAlignment {

    private:

        seq_len_t _pattern_start, _text_start, _len;
        AlignmentEvents _events;

    };


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