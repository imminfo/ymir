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


    /**
     *
     */
    struct AlignmentVectorBase {


        typedef std::vector<bool> events_storage_t;

            
        static const size_t default_start_reserve_size = 20;

        static const size_t default_events_reserve_size = 600;


        AlignmentVectorBase() {
            _starts.reserve(default_start_reserve_size);
            _events.reserve(default_events_reserve_size);
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


    /**
     *
     */
    struct NoGapAlignmentVector : public AlignmentVectorBase {

        NoGapAlignmentVector() : AlignmentVectorBase() {

        }


        /**
         * \brief Add a new alignment to the vector.
         */
        ///@{
        void addAlignment(seq_len_t p_start, seq_len_t t_start, seq_len_t size) {
            _data.push_back(p_start);
            _data.push_back(t_start);
            _data.push_back(size);
        }

        void addAlignment(seq_len_t p_start, seq_len_t t_start, const events_storage_t &vec) {
            _data.push_back(p_start);
            _data.push_back(t_start);
            _data.push_back(vec.size());
            _starts.push_back(_events.size());
            _events.insert(_events.end(), vec.begin(), vec.end());
        }
        ///@}


        bool isMismatch(seq_len_t i, seq_len_t j) const { 
#ifndef DNDEBUG
            if (_starts[i] + j >= _events.size()) {
                throw(std::runtime_error("Alignment vector: index is out of bounds."));
            }
#endif
            return _events[_starts[i] + j];
        }

    };


    /**
     *
     */
    struct GappedAlignmentVector : public AlignmentVectorBase {


        GappedAlignmentVector() : AlignmentVectorBase() {

        }


        void addAlignment(seq_len_t p_start, seq_len_t t_start, const events_storage_t &vec) {
            _data.push_back(p_start);
            _data.push_back(t_start);
            _data.push_back(vec.size());
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


    struct VDJAlignment {

        /**
         * \brief Move constructor for _segments, _alignments and _n_D_alignments.
         */
        VDJAlignment(NoGapAlignmentVector &&alignments) 
        {
        }


        VDJAlignment(const VDJAlignment &other) {

        }


        VDJAlignment& operator=(const VDJAlignment &other) {

        }


        virtual ~VDJAlignment()
        {
        }


        /** 
         * \brief Get the number of alignments for the specific gene.
         */
        ///@{
        seg_index_t nVar() const { return _segments[0]; }

        seg_index_t nJoi() const { return _segments[1]; }

        seg_index_t nDiv() const { return _segments[2]; }
        ///@}


        /**
         * \brief Get the index of the aligned gene segment for the specific gene.
         */
        ///@{
        seg_index_t getVar(size_t index) const { return _segments[3 + index]; }

        seg_index_t getJoi(size_t index) const { return _segments[3 + _segments[0] + index]; }

        seg_index_t getDiv(size_t index) const { return _segments[3 + _segments[0] + _segments[1] + index]; }
        ///@}


        /**
         * \brief Get alignments for the specific gene.
         */
        ///@{
        ///@}


        seq_len_t numDivAlignments(seg_index_t index) const {
            if (_n_D_alignments) {
                return _n_D_alignments[index];
            } else {
                return 0;
            }
        }

    protected:

        unique_ptr<seg_index_t[]> _segments;  /// Two concatenated vectors: vector of length 3 w/ numbers of aligned segments (V-J-D) and
                                              /// vector of indices of segments, aligned on this clone: V1--V2--V3--J1--J2--D1--D2--...
        NoGapAlignmentVector _alignments;  /// Vector of alignments for segments.
        unique_ptr<seq_len_t[]> _n_D_alignments;  /// Number of alignments for each D gene segment; vector's length == _segments[2]


        VDJAlignment();

    };

}

#endif