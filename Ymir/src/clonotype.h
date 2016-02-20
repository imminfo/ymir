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

#ifndef _CLONOTYPE_H
#define _CLONOTYPE_H


#include "alignment.h"


namespace ymir {
    

    struct Clonotype;


    /**
     * \typedef ClonotypePtr
     */
    typedef unique_ptr<Clonotype> ClonotypePtr;


    /**
    * \struct Clonotype
    */
    // struct Clonotype : public VDJAlignment {
    struct Clonotype {

        Clonotype(const std::string& sequence,
                  SequenceType seq_type,
                  Recombination recomb,
                  seg_index_t *segments,
                  seq_len_t *alignments,
                  seq_len_t *n_D_alignments = nullptr)
                : _sequence(sequence),
                  _seq_type(seq_type),
                  _recomb(recomb),
                  _segments(segments),
                  _alignments(alignments),
                  _n_D_alignments(n_D_alignments)
        { }


        Clonotype(const Clonotype& other)
                : _seq_type(other._seq_type),
                  _sequence(other._sequence),
                  _recomb(other._recomb)
        {
            _segments = new seg_index_t[3 + other._segments[0] + other._segments[1] + other._segments[2]];
            _segments[0] = other._segments[0];
            _segments[1] = other._segments[1];
            _segments[2] = other._segments[2];
            for (int i = 0; i < _segments[0] + _segments[1] + _segments[2]; ++i) {
                _segments[i + 3] = other._segments[i + 3];
            }

            int sum_align = 3*_segments[0] + 3*_segments[1];
            if (_segments[2]) {
                _n_D_alignments = new seq_len_t[_segments[2]];
                for  (int i = 0; i < _segments[2]; ++i) {
                    _n_D_alignments[i] = other._n_D_alignments[i];
                }
                for  (int i = 0; i < _segments[2]; ++i) {
                    sum_align += _n_D_alignments[i] * 3;
                }

            } else {
                _n_D_alignments = nullptr;
            }

            _alignments = new seq_len_t[sum_align];
            for (int i = 0; i < sum_align; ++i) {
                _alignments[i] = other._alignments[i];
            }
        }


        Clonotype(const ClonotypePtr& other)
                : _seq_type(other->_seq_type),
                  _sequence(other->_sequence),
                  _recomb(other->_recomb)
        {
            _segments = new seg_index_t[3 + other->_segments[0] + other->_segments[1] + other->_segments[2]];
            _segments[0] = other->_segments[0];
            _segments[1] = other->_segments[1];
            _segments[2] = other->_segments[2];
            for (int i = 0; i < _segments[0] + _segments[1] + _segments[2]; ++i) {
                _segments[i + 3] = other->_segments[i + 3];
            }

            int sum_align = 3*_segments[0] + 3*_segments[1];
            if (_segments[2]) {
                _n_D_alignments = new seq_len_t[_segments[2]];
                for  (int i = 0; i < _segments[2]; ++i) {
                    _n_D_alignments[i] = other->_n_D_alignments[i];
                }
                for  (int i = 0; i < _segments[2]; ++i) {
                    sum_align += _n_D_alignments[i] * 3;
                }

            } else {
                _n_D_alignments = nullptr;
            }

            _alignments = new seq_len_t[sum_align];
            for (int i = 0; i < sum_align; ++i) {
                _alignments[i] = other->_alignments[i];
            }
        }


        ///@{
        Clonotype& operator=(const Clonotype &other) {
            _seq_type = other._seq_type;
            _sequence = other._sequence;
            _recomb = other._recomb;

            _segments = new seg_index_t[3 + other._segments[0] + other._segments[1] + other._segments[2]];
            _segments[0] = other._segments[0];
            _segments[1] = other._segments[1];
            _segments[2] = other._segments[2];
            for (int i = 0; i < _segments[0] + _segments[1] + _segments[2]; ++i) {
                _segments[i + 3] = other._segments[i + 3];
            }

            int sum_align = 3*_segments[0] + 3*_segments[1];
            if (_segments[2]) {
                _n_D_alignments = new seq_len_t[_segments[2]];
                for  (int i = 0; i < _segments[2]; ++i) {
                    _n_D_alignments[i] = other._n_D_alignments[i];
                }
                for  (int i = 0; i < _segments[2]; ++i) {
                    sum_align += _n_D_alignments[i] * 3;
                }

            } else {
                _n_D_alignments = nullptr;
            }

            _alignments = new seq_len_t[sum_align];
            for (int i = 0; i < sum_align; ++i) {
                _alignments[i] = other._alignments[i];
            }
        }

        Clonotype& operator=(const ClonotypePtr &other) {
            _seq_type = other->_seq_type;

            _sequence = other->_sequence;
            _recomb = other->_recomb;

            _segments = new seg_index_t[3 + other->_segments[0] + other->_segments[1] + other->_segments[2]];
            _segments[0] = other->_segments[0];
            _segments[1] = other->_segments[1];
            _segments[2] = other->_segments[2];
            for (int i = 0; i < _segments[0] + _segments[1] + _segments[2]; ++i) {
                _segments[i + 3] = other->_segments[i + 3];
            }

            int sum_align = 3*_segments[0] + 3*_segments[1];
            if (_segments[2]) {
                _n_D_alignments = new seq_len_t[_segments[2]];
                for  (int i = 0; i < _segments[2]; ++i) {
                    _n_D_alignments[i] = other->_n_D_alignments[i];
                }
                for  (int i = 0; i < _segments[2]; ++i) {
                    sum_align += _n_D_alignments[i] * 3;
                }

            } else {
                _n_D_alignments = nullptr;
            }

            _alignments = new seq_len_t[sum_align];
            for (int i = 0; i < sum_align; ++i) {
                _alignments[i] = other->_alignments[i];
            }
        }
        ///@}


        virtual ~Clonotype() {
            if (_segments) {
                delete [] _segments;
            }
            if (_alignments) {
                delete [] _alignments;
            }
            if (_n_D_alignments) {
                delete [] _n_D_alignments;
            }
        }


        /**
         *
         */
        ///@{
        const sequence_t& sequence() const { return _sequence; }

        const sequence_t& nuc_sequence() const {
#ifndef DNDEBUG
            check_and_throw(_seq_type != NUCLEOTIDE, "Clonotype's call to nuc_sequence() is incorrect: wrong sequence type.");
#endif
            return _sequence;
        }

        const sequence_t& aa_sequence() const {
#ifndef DNDEBUG
            check_and_throw(_seq_type == UNDEF_SEQ_TYPE, "Clonotype's call to aa_sequence() is incorrect: undefined sequence type.");
#endif
            if (_seq_type == NUCLEOTIDE) {
                return translate(_sequence);
            } else {
                return _sequence;
            }
        }
        ///@}


        sequence_t::const_iterator seq_iterator(seq_len_t pos) const { return _sequence.cbegin() + pos; }


        Recombination recombination() const { return _recomb; }


        SequenceType sequence_type() const { return _seq_type; }


        /** Get index of aligned gene segment for specific gene. */
        ///@{
        seg_index_t getVar(size_t index) const { return _segments[3 + index]; }

        seg_index_t getJoi(size_t index) const { return _segments[3 + _segments[0] + index]; }

        seg_index_t getDiv(size_t index) const { return _segments[3 + _segments[0] + _segments[1] + index]; }
        ///@}


        /** Get number of alignments for specific gene. */
        ///@{
        seg_index_t nVar() const { return _segments[0]; }

        seg_index_t nJoi() const { return _segments[1]; }

        seg_index_t nDiv() const { return _segments[2]; }
        ///@}


        ///@{
        Alignment getVarAlignment(seg_index_t gene_index) const { return Alignment(_alignments + 3*gene_index); }

        Alignment getJoiAlignment(seg_index_t gene_index) const { return Alignment(_alignments + 3*_segments[0] + 3*gene_index); }

        Alignment getDivAlignment(seg_index_t gene_index, seg_index_t alignment_index) const {
            seq_len_t shift = 3*_segments[0] + 3*_segments[1];
            for (size_t i = 0; i < gene_index; ++i) {
                shift += 3 * _n_D_alignments[i];
            }
            return Alignment(_alignments + shift + 3*alignment_index);
        }
        ///@}


        seq_len_t numDivAlignments(seg_index_t index) const {
            if (_n_D_alignments) {
                return _n_D_alignments[index];
            } else {
                return 0;
            }
        }

    protected:

        Recombination _recomb;

        SequenceType _seq_type;

        std::string _sequence; //* CDR3 or full nucleotide or amino acid sequence of a clone. */

        seg_index_t *_segments; /// Two concatenated vectors: vector of length 3 w/ numbers of aligned segments (V-J-D) and
                                /// vector of indices of segments, aligned on this clone: V1--V2--V3--J1--J2--D1--D2--...

        seq_len_t *_alignments;  /// Vector of 3-tuples 1-based alignments (gene segment start, reference sequence start, alignment length)
                                 /// for alignments for each aligned gene segment.

        seq_len_t *_n_D_alignments; //* Number of aligned D gene segments (i.e., 3-tuples); vector's length == _segments[2] */


        // NoGapAlignmentVector _alignments;


        /**
         *
         */
        Clonotype() {}

    };

}

#endif