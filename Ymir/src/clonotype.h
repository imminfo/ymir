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


#include "types.h"

#include "vector"


namespace ymir {

    #define CLONOTYPEBUILDER_VSEG_DEFAULT_RESERVE_SIZE 3
    #define CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE 2
    #define CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE 40


    /**
    * \struct Clonotype
    */
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
                : _sequence(other._sequence),
                  _recomb(other._recomb)
        {
            _segments = new seg_index_t[3 + other._segments[0] + other._segments[1] + other._segments[2]];
            _segments[0] = other._segments[0];
            _segments[1] = other._segments[1];
            _segments[2] = other._segments[2];
            for (int i = 0; i < _segments[0] + _segments[1] + _segments[2]; ++i) {
                _segments[i + 3] = other._segments[i + 3];
            }

            int sum_align = _segments[0] + _segments[1];
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

        // Clonotype& operator=


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


        const std::string& sequence() const { return _sequence; }


        std::string::const_iterator seq_iterator(seq_len_t pos) const { return _sequence.cbegin() + pos; }


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


        seq_len_t getVend(seg_index_t index) const { return _alignments[index]; }


        seq_len_t getJstart(seg_index_t index) const { return _alignments[index + _segments[0]]; }


        Alignment getDivAlignment(seg_index_t Dseg_index, seg_index_t index) const {
            seq_len_t shift = _segments[0] + _segments[1];
            for (size_t i = 0; i < Dseg_index; ++i) {
                shift += 3 * _n_D_alignments[i];
            }
            return Alignment(_alignments + shift + 3*index);

            return this->getAlignment(Dseg_index, index, _segments[0] + _segments[1]);
        }

        seq_len_t numDivAlignments(seg_index_t index) const {
            if (_n_alignments) {
                return _n_alignments[index];
            } else {
                return 0;
            }
        }

    protected:

        Recombination _recomb;

        SequenceType _seq_type;

        std::string _sequence; //* CDR3 or full nucleotide or amino acid sequence of a clone. */

        seg_index_t *_segments; /// Three concatenated vectors: vector of length 3 w/ numbers of aligned segments (V-J-D) and
                                /// vector of indices of segments, aligned on this clone: V1--V2--V3--J1--J2--D1--D2--...

        seq_len_t *_alignments;  /// Vector of 3-tuples 1-based alignments (gene segment start, gene segment end, reference sequence start)
                                 /// for alignments for each aligned gene segment.

        seq_len_t *_n_alignments; //* Number of alignments (i.e., 3-tuples) for each aligned segment; vector's length == _segments[0] + _segments[1] + _segments[2] */


        /**
         *
         */
        Clonotype() {}


        Alignment getAlignment(seg_index_t gene_index,
                               seg_index_t alignment_index,
                               seg_index_t search_start,
                               seg_index_t search_end) {
            size_t shift;
            for (size_t i = 0; i < search_end; ++i) {
                shift += 3 * _n_alignments[i];
            }
            return Alignment(_alignments + shift + 3*alignment_index);
        }

    };


    /**
    * \class ClonotypeBuilder
    */
    class ClonotypeBuilder : protected Clonotype {
    public:

        ClonotypeBuilder() {
            reset();
        }


        virtual ~ClonotypeBuilder() {

        }


        /**
        * \brief Build clone alignment structure with stored information.
        *
        * \return Pointer to the newly created ClonotypeAlignment object.
        */
        Clonotype buildClonotype() {
//            segindex_t *_segments = new segindex_t[3];
            _segments = new seg_index_t[3 + _Vseg.size() + _Jseg.size() + _n_Dalign.size()];
            _segments[0] = _Vseg.size();
            _segments[1] = _Jseg.size();
            _segments[2] = _n_Dalign.size();

            for (int i = 0; i < _segments[0]; ++i) {
                _segments[3 + i] = _Vseg[i];
            }
            for (int i = _segments[0]; i < _segments[0] + _segments[1]; ++i) {
                _segments[3 + i] = _Jseg[i - _segments[0]];
            }
            int prev = 0;
            for (int i = _segments[0] + _segments[1]; i < _segments[0] + _segments[1] + _segments[2]; ++i) {
                _segments[3 + i] = _Dseg[prev];
                prev += _n_Dalign[i - _segments[0] - _segments[1]];
            }

            _alignments = new seq_len_t[_segments[0] + _segments[1] + _Dalign.size()];
            for (int i = 0; i < _segments[0]; ++i) {
                _alignments[i] = _Valign[i];
            }
            for (int i = _segments[0]; i < _segments[0] + _segments[1]; ++i) {
                _alignments[i] = _Jalign[i - _segments[0]];
            }

            if (_segments[2]) {
                for (int i = _segments[0] + _segments[1]; i < _segments[0] + _segments[1] + _Dalign.size(); ++i) {
                    _alignments[i] = _Dalign[i - _segments[0] - _segments[1]];
                }

                _n_D_alignments = new seq_len_t[_segments[2]];
                for (int i = 0; i < _n_Dalign.size(); ++i) {
                    _n_D_alignments[i] = _n_Dalign[i];
                }

            } else {
                _n_D_alignments = nullptr;
            }


            Clonotype cla = Clonotype(_sequence, _seq_type, _recomb, _segments, _alignments, _n_D_alignments);
            reset();
            return cla;
        }


        /**
         *
         */
        ///@{
        ClonotypeBuilder& setSequence(const std::string& seq) { this->_sequence = seq; return *this; }

        ClonotypeBuilder& setNucleotideSeq() { _seq_type = NUCLEOTIDE; return *this; }

        ClonotypeBuilder& setAminoAcidSeq() { _seq_type = AMINOACID; return *this; }

        ClonotypeBuilder& setRecombination(Recombination recomb) { _recomb = recomb; return *this; }
        ///@}


        ///@{
        ClonotypeBuilder& addVarAlignment(seg_index_t vseg, seq_len_t vstart, seq_len_t vend, seq_len_t seqstart) {
            return this->addAlignment(_Vseg, _Valign, _n_Valign, vseg, vstart, vend, seqstart);
        }

        ClonotypeBuilder& addJoiAlignment(seg_index_t jseg, seq_len_t jstart, seq_len_t jend, seq_len_t seqstart) {
            return this->addAlignment(_Jseg, _Jalign, _n_Jalign, jseg, jstart, jend, seqstart);
        }

        ClonotypeBuilder& addDivAlignment(seg_index_t dseg, seq_len_t dstart, seq_len_t dend, seq_len_t seqstart) {
            return this->addAlignment(_Dseg, _Dalign, _n_Dalign, dseg, dstart, dend, seqstart);
        }

        ClonotypeBuilder& addDivAlignment(seg_index_t dseg, const Alignment& dalignment) {
            return addDivAlignment(dseg, dalignment.gene_start(), dalignment.gene_end(), dalignment.seq_start());
        }
        ///@}


        /**
        * \brief Reset the builder to the initial state and remove all stored data for clone building.
        */
        void reset() {
            _segments = nullptr;
            _alignments = nullptr;
            _n_D_alignments = nullptr;

            _Vseg.clear();
            _Vseg.reserve(CLONOTYPEBUILDER_VSEG_DEFAULT_RESERVE_SIZE);
            _Valign.clear();
            _Valign.reserve(CLONOTYPEBUILDER_VSEG_DEFAULT_RESERVE_SIZE);

            _n_Valign.clear();

            _Jseg.clear();
            _Jseg.reserve(CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE);
            _Jalign.clear();
            _Jalign.reserve(CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE);

            _n_Jalign.clear();

            _Dseg.clear();
            _Dseg.reserve(CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE);
            _Dalign.clear();
            _Dalign.reserve(CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE * 4);

            _n_Dalign.clear();
        }

    protected:

        std::vector<seg_index_t> _Vseg, _Jseg, _Dseg; // gene segment aligned
        std::vector<seq_len_t> _Valign, _Jalign, _Dalign; // alignments for each gene segment
        std::vector<seg_index_t> _n_Valign, _n_Jalign, _n_Dalign; // number of D alignments


        ClonotypeBuilder& addAlignment(std::vector<seg_index_t> &seg_vec,
                                       std::vector<seq_len_t> &align_vec,
                                       std::vector<seg_index_t> n_align_vec,
                                       seg_index_t gene_index,
                                       seq_len_t gene_start,
                                       seq_len_t gene_end,
                                       seq_len_t seq_start)
        {
            // check for recombination type?

            if (seg_vec.size() == 0 || gene_index != seg_vec[seg_vec.size() - 1]) {
                n_align_vec.push_back(0);
            }

            n_align_vec[n_align_vec.size() - 1] += 1;

            seg_vec.push_back(gene_index);
            align_vec.push_back(gene_start);
            align_vec.push_back(gene_end);
            align_vec.push_back(seq_start);
            return *this;
        }

    };

}

#endif