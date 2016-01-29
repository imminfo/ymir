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

    #define CLONOTYPEBUILDER_VSEG_DEFAULT_RESERVE_SIZE 9
    #define CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE 6
    #define CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE 30


    struct Clonotype;


    /**
     * \typedef ClonotypePtr
     */
    typedef unique_ptr<Clonotype> ClonotypePtr;


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
                : _sequence(other->_sequence),
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


        /**
         *
         */
        Clonotype() {}

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
        ClonotypePtr buildClonotype() {
            _segments = new seg_index_t[3 + _Vseg.size() + _Jseg.size() + _n_Dalign.size()];
            _segments[0] = _Vseg.size();
            _segments[1] = _Jseg.size();
            _segments[2] = _n_Dalign.size();

            for (int i = 0; i < _segments[0]; ++i) {
                _segments[3 + i] = _Vseg[i];
            }

            for (int i = 0; i < _segments[1]; ++i) {
                _segments[3 + i + _segments[0]] = _Jseg[i];
            }

            int prev = 0;
            for (int i = _segments[0] + _segments[1]; i < _segments[0] + _segments[1] + _segments[2]; ++i) {
                _segments[3 + i] = _Dseg[prev];
                prev += _n_Dalign[i - _segments[0] - _segments[1]];
            }

            _alignments = new seq_len_t[_Valign.size() + _Jalign.size() + _Dalign.size()];
            for (int i = 0; i < _Valign.size(); ++i) {
                _alignments[i] = _Valign[i];
            }

            for (int i = 0; i < _Jalign.size(); ++i) {
                _alignments[i + _Valign.size()] = _Jalign[i];
            }

            if (_segments[2]) {
                for (int i = 0; i < _Dalign.size(); ++i) {
                    _alignments[i + _Valign.size() + _Jalign.size()] = _Dalign[i];
                }

                _n_D_alignments = new seq_len_t[_segments[2]];
                for (int i = 0; i < _n_Dalign.size(); ++i) {
                    _n_D_alignments[i] = _n_Dalign[i];
                }

            } else {
                _n_D_alignments = nullptr;
            }

            ClonotypePtr cla(new Clonotype(_sequence, _seq_type, _recomb, _segments, _alignments, _n_D_alignments));
            reset();
            return std::move(cla);
        }


        /**
         *
         */
        ///@{
        ClonotypeBuilder& setSequence(const std::string& seq) { this->_sequence = seq; return *this; }

        ClonotypeBuilder& setSequenceType(SequenceType seq_type) { _seq_type = seq_type; return *this; }

        ClonotypeBuilder& setNucleotideSeq() { _seq_type = NUCLEOTIDE; return *this; }

        ClonotypeBuilder& setAminoAcidSeq() { _seq_type = AMINOACID; return *this; }

        ClonotypeBuilder& setRecombination(Recombination recomb) { _recomb = recomb; return *this; }
        ///@}


        /**
         * \brief Interface to add alignments to the builder.
         */
        ///@{
        ClonotypeBuilder& addVarAlignment(seg_index_t vseg, seq_len_t vstart, seq_len_t seqstart, seq_len_t alignment_len) {
            _Vseg.push_back(vseg);
            _Valign.push_back(vstart);
            _Valign.push_back(seqstart);
            _Valign.push_back(alignment_len);
            return *this;
        }

        ClonotypeBuilder& addJoiAlignment(seg_index_t jseg, seq_len_t jstart, seq_len_t seqstart, seq_len_t alignment_len) {
            _Jseg.push_back(jseg);
            _Jalign.push_back(jstart);
            _Jalign.push_back(seqstart);
            _Jalign.push_back(alignment_len);
            return *this;
        }

        ClonotypeBuilder& addDivAlignment(seg_index_t dseg, seq_len_t dstart, seq_len_t seqstart, seq_len_t alignment_len) {
            // check for recombination type?

            if (_Dseg.size() == 0 || dseg != _Dseg[_Dseg.size() - 1]) {
                _n_Dalign.push_back(0);
            }

            _n_Dalign[_n_Dalign.size() - 1] += 1;

            _Dseg.push_back(dseg);
            _Dalign.push_back(dstart);
            _Dalign.push_back(seqstart);
            _Dalign.push_back(alignment_len);
            return *this;
        }

        ClonotypeBuilder& addDivAlignment(seg_index_t dseg, const Alignment& dalignment) {
            return addDivAlignment(dseg, dalignment.gene_start(), dalignment.seq_start(), dalignment.length());
        }


        ClonotypeBuilder& addAlignment(GeneSegments gene_segment, seg_index_t seg, seq_len_t gene_start, seq_len_t seq_start, seq_len_t alignment_len) {
            switch (gene_segment) {

                case VARIABLE: {
                    this->addVarAlignment(seg, gene_start, seq_start, alignment_len);
                    break;
                }

                case DIVERSITY: {
                    this->addDivAlignment(seg, gene_start, seq_start, alignment_len);
                    break;
                }

                case JOINING: {
                    this->addJoiAlignment(seg, gene_start, seq_start, alignment_len);
                    break;
                }

                default: {}
            }
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

            _Jseg.clear();
            _Jseg.reserve(CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE);
            _Jalign.clear();
            _Jalign.reserve(CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE);

            _Dseg.clear();
            _Dseg.reserve(CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE);
            _Dalign.clear();
            _Dalign.reserve(CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE * 4);

            _n_Dalign.clear();
        }

    protected:

        // vector -> map so can insert segments in random order
        std::vector<seg_index_t> _Vseg, _Jseg, _Dseg; // gene segment aligned
        std::vector<seq_len_t> _Valign, _Jalign, _Dalign; // alignments for each gene segment
        std::vector<seg_index_t> _n_Dalign; // number of D alignments

    };

}

#endif