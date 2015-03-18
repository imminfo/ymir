/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vadim dot nazarov at mailbox dot com>
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

        Clonotype(const string& sequence, bool nucleotide, segindex_t *segments, segindex_t *n_segments,
                seq_len_t *alignments, seq_len_t *n_D_alignments = nullptr) : seq(sequence), nucleotide(nucleotide) {

            this->segments = segments;
            this->n_segments = n_segments;
            this->alignments = alignments;
            this->n_D_alignments = n_D_alignments;
        }


        Clonotype(const Clonotype& other) {
            this->seq = other.seq;
            this->nucleotide = other.nucleotide;

            this->n_segments = new segindex_t[3];
            this->n_segments[0] = other.n_segments[0];
            this->n_segments[1] = other.n_segments[1];
            this->n_segments[2] = other.n_segments[2];

            this->segments = new segindex_t[this->n_segments[0] + this->n_segments[1] + this->n_segments[2]];
            for (int i = 0; i < this->n_segments[0] + this->n_segments[1] + this->n_segments[2]; ++i) {
                this->segments[i] = other.segments[i];
            }

            int sum_align = this->n_segments[0] + this->n_segments[1];
            if (this->n_segments[2]) {
                this->n_D_alignments = new seq_len_t[this->n_segments[2]];
                for  (int i = 0; i < this->n_segments[2]; ++i) {
                    this->n_D_alignments[i] += other.n_D_alignments[i];
                }
                for  (int i = 0; i < this->n_segments[2]; ++i) {
                    sum_align += this->n_D_alignments[i] * 4;
                }
                for (int i = 0; i < sum_align; ++i) {
                    this->alignments[i] = other.alignments[i];
                }
            } else {
                this->n_D_alignments = nullptr;
            }
        }


        virtual ~Clonotype() {
            if (this->segments) {
                delete [] this->segments;
            }
            if (this->n_segments) {
                delete [] this->n_segments;
            }
            if (this->alignments) {
                delete [] this->alignments;
            }
            if (this->n_D_alignments) {
                delete [] this->n_D_alignments;
            }
        }


        const string& sequence() const { return this->seq; }


        bool is_nucleotide() const { return this->nucleotide; }


        /**
        * \brief Check if this clone was assembled in VDJ-recombination, i.e., it has D(iversity) gene segment.
        *
        * \return True if clone was assembled in VDJ-recombination; false otherwise.
        */
        bool is_vdj() const { return this->n_segments[2] > 0; }


        ///@{
        /** Get index of aligned gene segment for specific gene. */
        segindex_t getV(size_t index) const { return this->segments[index]; }
        segindex_t getJ(size_t index) const { return this->segments[this->n_segments[0] + index]; }
        segindex_t getD(size_t index) const { return this->segments[this->n_segments[0] + this->n_segments[1] + index]; }
        ///@}


        ///@{
        /** Get number of alignments for specific gene. */
        seq_len_t nV() const { return this->n_segments[0]; }
        seq_len_t nJ() const { return this->n_segments[1]; }
        seq_len_t nD() const { return this->n_segments[2]; }
        ///@}


        seq_len_t getVend(segindex_t index) const { return this->alignments[index]; }


        seq_len_t getJstart(segindex_t index) const { return this->alignments[index + this->n_segments[0]]; }


        d_alignment_t getDalignment(segindex_t Dseg_index, segindex_t index) const {
            seq_len_t shift = this->n_segments[0] + this->n_segments[1];
            for (size_t i = 0; i < Dseg_index; ++i) {
                shift += 4 * this->n_D_alignments[i];
            }
            return d_alignment_t(this->alignments + shift + 4*index);
        }

        seq_len_t nDalignments(segindex_t index) const {
            if (this->n_D_alignments) {
                return this->n_D_alignments[index];
            } else {
                return 0;
            }
        }

    protected:

        string seq; //* CDR3 or full nucleotide or amino acid sequence of a clone. */
        bool nucleotide;
        segindex_t *segments; //* Vector of indices of segments, aligned on this clone: V1--V2--V3--J1--J2--D1--D2--... */
        segindex_t *n_segments; //* Vector of length 3 w/ numbers aligned segments. 1st element - number of aligned V segments, and so on. */
        seq_len_t *alignments; //* Vector of alignments for the clone. 1st N elements is V ends, than M elements is J starts, each other are 4-tuples for Ds alignment - (D start, D end, seq start, seq end) */
        seq_len_t *n_D_alignments; //* Number of alignments (i.e., 4-tuples) for each aligned D segment; vector's length == n_segments[2] */


        Clonotype() {}

    };


    /**
    * \class ClonotypeBuilder
    */
    class ClonotypeBuilder : protected Clonotype {
    public:

        ClonotypeBuilder() {
            this->reset();
        }


        virtual ~ClonotypeBuilder() {
        }


        /**
        * \brief Build clone alignment structure with stored information.
        *
        * \return Pointer to the newly created ClonotypeAlignment object.
        */
        Clonotype buildClonotype() {
//            segindex_t *n_segments = new segindex_t[3];
            n_segments = new segindex_t[3];
            n_segments[0] = this->_Vseg.size();
            n_segments[1] = this->_Jseg.size();
            n_segments[2] = this->_n_Dalign.size();

//            segindex_t *segments = new segindex_t[n_segments[0] + n_segments[1] + n_segments[2]];
            segments = new segindex_t[n_segments[0] + n_segments[1] + n_segments[2]];
            for (int i = 0; i < n_segments[0]; ++i) {
                segments[i] = this->_Vseg[i];
            }
            for (int i = n_segments[0]; i < n_segments[0] + n_segments[1]; ++i) {
                segments[i] = this->_Jseg[i - n_segments[0]];
            }
            int prev = 0;
            for (int i = n_segments[0] + n_segments[1]; i < n_segments[0] + n_segments[1] + n_segments[2]; ++i) {
                segments[i] = this->_Dseg[i - n_segments[0] - n_segments[1] + prev];
                prev += this->_n_Dalign[i - n_segments[0] - n_segments[1]];
            }

            alignments = new seq_len_t[n_segments[0] + n_segments[1] + this->_Dalign.size()];
            for (int i = 0; i < n_segments[0]; ++i) {
                alignments[i] = this->_Valign[i];
            }
            for (int i = n_segments[0]; i < n_segments[0] + n_segments[1]; ++i) {
                alignments[i] = _Jalign[i - n_segments[0]];
            }

            if (n_segments[2]) {
                for (int i = n_segments[0] + n_segments[1]; i < n_segments[0] + n_segments[1] + _Dalign.size(); ++i) {
                    alignments[i] = _Dalign[i - n_segments[0] - n_segments[1]];
                }

                n_D_alignments = new seq_len_t[n_segments[2]];
                for (int i = 0; i < _n_Dalign.size(); ++i) {
                    n_D_alignments[i] = _n_Dalign[i];
                }
            } else {
                n_D_alignments = nullptr;
            }


            Clonotype cla = Clonotype(this->seq, this->nucleotide, segments, n_segments, alignments, n_D_alignments);
//            n_segments = nullptr;
//            segments = nullptr;
//            alignments = nullptr;
//            n_D_alignments = nullptr;
            this->reset();
            return cla;
        }


        /**
        * \brief Set sequence for building a clone.
        *
        * \param seq Clonotype sequence.
        */
        void setSequence(const string& seq) { this->seq = seq; }


        ///@{
        void setNucleotideSeq() { this->nucleotide = true; }
        void setAminoAcidSeq() { this->nucleotide = false; }
        ///@}


        ///@{
        void addValignment(segindex_t vseg, seq_len_t vend) {
            this->_Vseg.push_back(vseg);
            this->_Valign.push_back(vend);
        }
        void addJalignment(segindex_t jseg, seq_len_t jstart) {
            this->_Jseg.push_back(jseg);
            this->_Jalign.push_back(jstart);
        }
        void addDalignment(segindex_t dseg, seq_len_t dstart, seq_len_t dend, seq_len_t seqstart, seq_len_t seqend) {
            if (this->_Dseg.size() == 0 || dseg != this->_Dseg[this->_Dseg.size() - 1]) {
                this->_n_Dalign.push_back(0);
            }

            this->_n_Dalign[this->_n_Dalign.size() - 1] += 1;

            this->_Dseg.push_back(dseg);
            this->_Dalign.push_back(dstart);
            this->_Dalign.push_back(dend);
            this->_Dalign.push_back(seqstart);
            this->_Dalign.push_back(seqend);
        }
        void addDalignment(segindex_t dseg, const d_alignment_t& dalignment) {
            this->addDalignment(dseg, dalignment.Dstart, dalignment.Dend, dalignment.seqstart, dalignment.seqend);
        }
        ///@}


        /**
        * \brief Reset the builder to the initial state and remove all stored data for clone building.
        */
        void reset() {
            n_segments = nullptr;
            segments = nullptr;
            alignments = nullptr;
            n_D_alignments = nullptr;

            this->_Vseg.clear();
            this->_Vseg.reserve(CLONOTYPEBUILDER_VSEG_DEFAULT_RESERVE_SIZE);
            this->_Valign.clear();
            this->_Valign.reserve(CLONOTYPEBUILDER_VSEG_DEFAULT_RESERVE_SIZE);

            this->_Jseg.clear();
            this->_Jseg.reserve(CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE);
            this->_Jalign.clear();
            this->_Jalign.reserve(CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE);

            this->_Dseg.clear();
            this->_Dseg.reserve(CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE);
            this->_Dalign.clear();
            this->_Dalign.reserve(CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE * 4);
            this->_n_Dalign.clear();
        }

    protected:

        vector<segindex_t> _Vseg, _Jseg, _Dseg; // gene segment aligned
        vector<seq_len_t> _Valign, _Jalign, _Dalign; // alignments for each gene segment
        vector<segindex_t> _n_Dalign; // number of D alignments

    };


    /**
    * \struct ClonotypeMetadata
    *
    * \brief Metadata clone: indices of events corresponding to
    * gene segments' alignment, deletions and nucleotide insertions.
    */
    struct ClonotypeMetadata {

    public:

        ClonotypeMetadata(const Clonotype& clone) {

        }


        ~ClonotypeMetadata() {

        }

    private:

        bool _nucleotide, _vdj;

        string _sequence;

        eventind_t *_v_indices, *_j_indices, *_d_indices; //* Indices of event families of aligned gene segments (family event indices). */
        seq_len_t *_v_dels, *_j_dels, *_d_dels; //* Number of deletions for each aligned gene segments. */

        seq_len_t *_n_segments, *_n_D_segments; //* Number of aligned gene segments. */

        eventind_t *_ins_lens; //* Indices of families of insertions lengths. Size = _n_segments - 1 */

    };
}

#endif