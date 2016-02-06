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

#ifndef _CLONOTYPE_BUILDER_H
#define _CLONOTYPE_BUILDER_H


#include "clonotype.h"
#include "vdj_alignment_builder.h"


namespace ymir {

    #define CLONOTYPEBUILDER_VSEG_DEFAULT_RESERVE_SIZE 9
    #define CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE 6
    #define CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE 30


    /**
    * \class ClonotypeBuilder
    */
    // class ClonotypeBuilder : protected Clonotype, public VDJAlignmentBuilder {
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
            _seq_type = UNDEF_SEQ_TYPE;

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