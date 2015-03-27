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


using namespace std;


namespace ymir {

    #define CLONOTYPEBUILDER_VSEG_DEFAULT_RESERVE_SIZE 3
    #define CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE 2
    #define CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE 40


    /**
    * \struct Clonotype
    */
    struct Clonotype {

        Clonotype(const string& sequence, bool nucleotide, segindex_t *segments,
                seq_len_t *alignments, seq_len_t *n_D_alignments = nullptr) : seq(sequence), nucleotide(nucleotide) {

            this->segments = segments;
            this->alignments = alignments;
            this->n_D_alignments = n_D_alignments;
        }


        Clonotype(const Clonotype& other) {
            seq = other.seq;
            nucleotide = other.nucleotide;

            segments = new segindex_t[3 + other.segments[0] + other.segments[1] + other.segments[2]];
            segments[0] = other.segments[0];
            segments[1] = other.segments[1];
            segments[2] = other.segments[2];
            for (int i = 0; i < segments[0] + segments[1] + segments[2]; ++i) {
                segments[i + 3] = other.segments[i + 3];
            }

            int sum_align = segments[0] + segments[1];
            if (segments[2]) {
                n_D_alignments = new seq_len_t[segments[2]];
                for  (int i = 0; i < segments[2]; ++i) {
                    n_D_alignments[i] += other.n_D_alignments[i];
                }
                for  (int i = 0; i < segments[2]; ++i) {
                    sum_align += n_D_alignments[i] * 4;
                }
            } else {
                n_D_alignments = nullptr;
            }

            alignments = new seq_len_t[sum_align];
            for (int i = 0; i < sum_align; ++i) {
                alignments[i] = other.alignments[i];
            }
        }


        virtual ~Clonotype() {
            if (segments) {
                delete [] segments;
            }
            if (alignments) {
                delete [] alignments;
            }
            if (n_D_alignments) {
                delete [] n_D_alignments;
            }
        }


        const string& sequence() const { return seq; }

        string::iterator seq_iterator(seq_len_t pos) { return seq.begin() + pos; }


        bool is_nucleotide() const { return nucleotide; }


        /**
        * \brief Check if this clone was assembled in VDJ-recombination, i.e., it has D(iversity) gene segment.
        *
        * \return True if clone was assembled in VDJ-recombination; false otherwise.
        */
        bool is_vdj() const { return segments[2] > 0; }


        ///@{
        /** Get index of aligned gene segment for specific gene. */
        segindex_t getVar(size_t index) const { return segments[3 + index]; }
        segindex_t getJoi(size_t index) const { return segments[3 + segments[0] + index]; }
        segindex_t getDiv(size_t index) const { return segments[3 + segments[0] + segments[1] + index]; }
        ///@}


        ///@{
        /** Get number of alignments for specific gene. */
        segindex_t nVar() const { return segments[0]; }
        segindex_t nJoi() const { return segments[1]; }
        segindex_t nDiv() const { return segments[2]; }
        ///@}


        seq_len_t getVend(segindex_t index) const { return alignments[index]; }


        seq_len_t getJstart(segindex_t index) const { return alignments[index + segments[0]]; }


        d_alignment_t getDalignment(segindex_t Dseg_index, segindex_t index) const {
            seq_len_t shift = segments[0] + segments[1];
            for (size_t i = 0; i < Dseg_index; ++i) {
                shift += 4 * n_D_alignments[i];
            }
            return d_alignment_t(alignments + shift + 4*index);
        }

        seq_len_t nDalignments(segindex_t index) const {
            if (n_D_alignments) {
                return n_D_alignments[index];
            } else {
                return 0;
            }
        }

    protected:

        string seq; //* CDR3 or full nucleotide or amino acid sequence of a clone. */
        bool nucleotide;
        segindex_t *segments; /// Two concatenated vectors: vector of length 3 w/ numbers of aligned segments (V-J-D) and
                              /// vector of indices of segments, aligned on this clone: V1--V2--V3--J1--J2--D1--D2--...
        seq_len_t *alignments; //* Vector of 1-based alignments for the clone. 1st N elements is V ends, than M elements is J starts, each other are 4-tuples for Ds alignment - (D start, D end, seq start, seq end) */
        seq_len_t *n_D_alignments; //* Number of alignments (i.e., 4-tuples) for each aligned D segment; vector's length == segments[2] */


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
        Clonotype buildClonotype() {
//            segindex_t *segments = new segindex_t[3];
            segments = new segindex_t[3 + _Vseg.size() + _Jseg.size() + _n_Dalign.size()];
            segments[0] = _Vseg.size();
            segments[1] = _Jseg.size();
            segments[2] = _n_Dalign.size();

            for (int i = 0; i < segments[0]; ++i) {
                segments[3 + i] = _Vseg[i];
            }
            for (int i = segments[0]; i < segments[0] + segments[1]; ++i) {
                segments[3 + i] = _Jseg[i - segments[0]];
            }
            int prev = 0;
            for (int i = segments[0] + segments[1]; i < segments[0] + segments[1] + segments[2]; ++i) {
                segments[3 + i] = _Dseg[i - segments[0] - segments[1] + prev];
                prev += _n_Dalign[i - segments[0] - segments[1]];
            }

            alignments = new seq_len_t[segments[0] + segments[1] + _Dalign.size()];
            for (int i = 0; i < segments[0]; ++i) {
                alignments[i] = _Valign[i];
            }
            for (int i = segments[0]; i < segments[0] + segments[1]; ++i) {
                alignments[i] = _Jalign[i - segments[0]];
            }

            if (segments[2]) {
                for (int i = segments[0] + segments[1]; i < segments[0] + segments[1] + _Dalign.size(); ++i) {
                    alignments[i] = _Dalign[i - segments[0] - segments[1]];
                }

                n_D_alignments = new seq_len_t[segments[2]];
                for (int i = 0; i < _n_Dalign.size(); ++i) {
                    n_D_alignments[i] = _n_Dalign[i];
                }
            } else {
                n_D_alignments = nullptr;
            }


            Clonotype cla = Clonotype(seq, nucleotide, segments, alignments, n_D_alignments);
            reset();
            return cla;

//            return value optimisation:
//            return Clonotype(seq, nucleotide, segments, alignments, n_D_alignments);
        }


        /**
        * \brief Set sequence for building a clone.
        *
        * \param seq Clonotype sequence.
        */
        ClonotypeBuilder& setSequence(const string& seq) { this->seq = seq; return *this; }


        ///@{
        ClonotypeBuilder& setNucleotideSeq() { nucleotide = true; return *this; }
        ClonotypeBuilder& setAminoAcidSeq() { nucleotide = false; return *this; }
        ///@}


        ///@{
        ClonotypeBuilder& addValignment(segindex_t vseg, seq_len_t vend) {
            _Vseg.push_back(vseg);
            _Valign.push_back(vend);
            return *this;
        }
        ClonotypeBuilder& addJalignment(segindex_t jseg, seq_len_t jstart) {
            _Jseg.push_back(jseg);
            _Jalign.push_back(jstart);
            return *this;
        }
        ClonotypeBuilder& addDalignment(segindex_t dseg, seq_len_t dstart, seq_len_t dend, seq_len_t seqstart, seq_len_t seqend) {
            if (_Dseg.size() == 0 || dseg != _Dseg[_Dseg.size() - 1]) {
                _n_Dalign.push_back(0);
            }

            _n_Dalign[_n_Dalign.size() - 1] += 1;

            _Dseg.push_back(dseg);
            _Dalign.push_back(dstart);
            _Dalign.push_back(dend);
            _Dalign.push_back(seqstart);
            _Dalign.push_back(seqend);
            return *this;
        }
        ClonotypeBuilder& addDalignment(segindex_t dseg, const d_alignment_t& dalignment) {
            addDalignment(dseg, dalignment.Dstart, dalignment.Dend, dalignment.seqstart, dalignment.seqend);
            return *this;
        }
        ///@}


        /**
        * \brief Reset the builder to the initial state and remove all stored data for clone building.
        */
        void reset() {
            segments = nullptr;
            alignments = nullptr;
            n_D_alignments = nullptr;

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

        vector<segindex_t> _Vseg, _Jseg, _Dseg; // gene segment aligned
        vector<seq_len_t> _Valign, _Jalign, _Dalign; // alignments for each gene segment
        vector<segindex_t> _n_Dalign; // number of D alignments

    };

}

#endif