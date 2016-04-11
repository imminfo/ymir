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

#ifndef _ALIGNER_H_
#define _ALIGNER_H_


#include "aligner_parameters.h"
#include "alignment_matrix.h"
#include "clonotype_builder.h"
#include "genesegment.h"


using namespace std;


namespace ymir {

    /**
     * \class VDJAlignerBase
     */
    template <typename ClonotypeType>
    class VDJAlignerBase : public ClonotypeBuilder<ClonotypeType> {

    public:


        typedef ClonotypeType clonotype_t;


        typedef typename ClonotypeType::vdj_alignment_t vdj_alignment_t;
        
        
        typedef typename vdj_alignment_t::alignment_vector_t alignment_vector_t;


        /**
         *
         */
        VDJAlignerBase()
        {
        }


        /**
         *
         */
        VDJAlignerBase(const VDJRecombinationGenes &genes,
                       VDJAlignerParameters params)
                : _genes(genes), _params(params)
        {
        }


        virtual ~VDJAlignerBase()
        {
        }


        void set_genes(const VDJRecombinationGenes &genes) { _genes = genes; }


        void set_parameters(VDJAlignerParameters params) { _params = params; }


        /**
         * \brief General methods for alignment.
         */
        ///@{
        alignment_vector_t alignVar(seg_index_t id, const sequence_t &sequence) const {
            alignment_vector_t vec;
            this->_alignVar(id, _genes.V()[id].sequence, _sequence, &vec);
            return vec;
        }

        alignment_vector_t alignDiv(seg_index_t id, const sequence_t &sequence) const {
            alignment_vector_t vec;
            this->_alignDiv(id, _genes.D()[id].sequence, _sequence, &vec);
            return vec;
        }

        alignment_vector_t alignJoi(seg_index_t id, const sequence_t &sequence) const {
            alignment_vector_t vec;
            this->_alignJoi(id, _genes.J()[id].sequence, _sequence, &vec);
            return vec;
        }
        ///@}


        /**
         * \brief Align the given sequence to all gene segments of the specific gene.
         *
         * \param sequence sequence sequence.
         *
         * \return True if has been aligned at least one gene segment, False if no gene segments have been aligned.
         */
        ///@{
        bool alignVar() {
            alignment_vector_t vec;
            for (seg_index_t id = 1; id <= _genes.V().max(); ++id) {
                this->_alignVar(id, _genes.V()[id].sequence, _sequence, &vec);
            }
            this->addVarAlignment(vec);
            return vec.size() != 0;
        }

        bool alignDiv() {
            alignment_vector_t vec;
            for (seg_index_t id = 1; id <= _genes.D().max(); ++id) {
                vec.clear(); // TODO: some strange behaviour here after I added this line. Be careful. Maybe there is some bug in clear().
                this->_alignDiv(id, _genes.D()[id].sequence, _sequence, &vec);
                if (vec.size()) {
                    this->addDivAlignment(vec);
                }
            }
            return vec.size() != 0;
        }

        bool alignJoi() {
            alignment_vector_t vec;
            for (seg_index_t id = 1; id <= _genes.J().max(); ++id) {
                this->_alignJoi(id, _genes.J()[id].sequence, _sequence, &vec);
            }
            this->addJoiAlignment(vec);
            return vec.size() != 0;
        }
        ///@}


    protected:

        VDJAlignerParameters _params;
        VDJRecombinationGenes _genes;


        /**
         * \brief Internal alignment functions. Inherit from this functions
         * to perform alignment.
         *
         * \param id
         * \param gene
         * \param sequence
         * \param vec
         */
        ///@{
        virtual void _alignVar(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, alignment_vector_t *avec) const = 0;

        virtual void _alignDiv(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, alignment_vector_t *avec) const = 0;

        virtual void _alignJoi(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, alignment_vector_t *avec) const = 0;
        ///@}

    };


    /**
     * \class NaiveCDR3NucleotideAligner
     *
     * \brief CDR3-only alignment without errors - align V starting from the left edge,
     * J starting from the right edge, and align D everywhere.
     */
    class NaiveCDR3NucleotideAligner : public VDJAlignerBase<ClonotypeNuc> {
    public:

        NaiveCDR3NucleotideAligner()
        {
        }


        NaiveCDR3NucleotideAligner(const VDJRecombinationGenes &genes,
                                   VDJAlignerParameters params)
            : VDJAlignerBase<ClonotypeNuc>(genes, params)
        {
        }


    protected:

        /**
         *
         */
        ///@{
        void _alignVar(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0;

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[i] != text[i]) { break; }
                matches += 1;
            }

            avec->addAlignment(gene, 1, 1, matches);
        }

        void _alignDiv(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t match_min_len = _params.min_D_len;
            seq_len_t t_size = text.size(), p_size = pattern.size(), min_size = min(t_size, p_size), min_subsize;
            bool open_match;
            seq_len_t p_start, t_start;

            for (seq_len_t pattern_i = 0 /* WTF?! */ ; pattern_i < p_size - match_min_len + 1; ++pattern_i) {
                open_match = false;
                min_subsize = min(p_size - pattern_i, (int) t_size);
                for (seq_len_t i = 0; i < min_subsize; ++i) {
                    if (pattern[pattern_i + i] == text[i]) {
                        if (!open_match) {
                            p_start = pattern_i + i;
                            t_start = i;
                            open_match = true;
                        }
                    } else if (open_match) {
                        if ((pattern_i + i - p_start) >= match_min_len) {
                            avec->addAlignment(gene, p_start + 1, t_start + 1, pattern_i + i - p_start);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (pattern_i + min_subsize - p_start) >= match_min_len) {
                    avec->addAlignment(gene, p_start + 1, t_start + 1, pattern_i + min_subsize - p_start);
                }
            }

            for (seq_len_t text_i = 1; text_i < t_size - match_min_len + 1; ++text_i) {
                open_match = false;
                min_subsize = min((int) p_size, t_size - text_i);
                for (seq_len_t i = 0; i < min_subsize; ++i) {
                    if (pattern[i] == text[text_i + i]) {
                        if (!open_match) {
                            p_start = i;
                            t_start = text_i + i;
                            open_match = true;
                        }
                    } else if (open_match) {
                        if ((i - p_start) >= match_min_len) {
                            avec->addAlignment(gene, p_start + 1, t_start + 1, i - p_start);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (min_subsize - p_start) >= match_min_len) {
                    avec->addAlignment(gene, p_start + 1, t_start + 1, min_subsize - p_start);
                }
            }
        }

        void _alignJoi(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0;

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[p_size - i - 1] != text[t_size - i - 1]) { break; }
                matches += 1;
            }

            avec->addAlignment(gene, p_size - matches + 1, t_size - matches + 1, matches);
        }
        ///@}

    };


    //
    // CDR3-only alignment with errors - align V starting from the left edge, 
    // J starting from the right edge, and align D everywhere.
    //
    class CDR3NucleotideAligner : public VDJAlignerBase<ClonotypeNuc> {
    public:

        CDR3NucleotideAligner()
        {
        }


        CDR3NucleotideAligner(const VDJRecombinationGenes &genes,
                              VDJAlignerParameters _params)
                : VDJAlignerBase<ClonotypeNuc>(genes, _params)
        {
        }


    protected:

        /**
         *
         */
        ///@{
        void _alignVar(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size();
            NoGapAlignment::events_storage_t vec;
            vec.reserve(min(p_size, t_size) + 1);
            alignment_score_t score = 0, val; //, max_score = 0;

            vec.push_back(pattern[0] != text[0]);
            score += pattern[0] == text[0] ? _params.score.v_score.match : _params.score.v_score.mism;
            for (seq_len_t i = 1; i < min(p_size, t_size); ++i) {
                vec.push_back(pattern[i] != text[i]);
//                val = pattern[i] == text[i] ? _params.score.v_score.match : (_params.score.v_score.mism - _params.score.v_score.acc_mism*(pattern[i] != text[i]));
                score += pattern[i] == text[i] ? (_params.score.v_score.match + _params.score.v_score.acc_match * (pattern[i - 1] == text[i - 1])) : _params.score.v_score.mism;
                // max_score = std::max(max_score, score);
            }

            // if (max_score >= _params.threshold.v_threshold) {
            if (score >= _params.threshold.v_threshold) {
                avec->addAlignment(gene, 1, 1, vec);
            }
        }

        void _alignDiv(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t match_min_len = _params.min_D_len;
            seq_len_t t_size = text.size(), p_size = pattern.size(), min_size = min(t_size, p_size), min_subsize;
            seq_len_t p_start, t_start;
            AlignmentVectorBase::events_storage_t bitvec;
            bitvec.reserve(p_size + 1);

            for (seq_len_t pattern_i = 0; pattern_i < p_size - match_min_len + 1; ++pattern_i) {
                min_subsize = min(p_size - pattern_i, (int) t_size);
                if (min_subsize >= match_min_len) {
                    bitvec.resize(min_subsize);
                    for (seq_len_t i = 0; i < min_subsize; ++i) {
                        bitvec[i] = pattern[pattern_i + i] != text[i];
                    }
                    avec->addAlignment(gene, pattern_i + 1, 1, bitvec);
                }
            }

            for (seq_len_t text_i = 1; text_i < t_size - match_min_len + 1; ++text_i) {
                min_subsize = min((int) p_size, t_size - text_i);
                if (min_subsize >= match_min_len) {
                    bitvec.resize(min_subsize);
                    for (seq_len_t i = 0; i < min_subsize; ++i) {
                        bitvec[i] = pattern[i] != text[text_i + i];
                    }
                    avec->addAlignment(gene, 1, text_i + 1, bitvec);
                }
            }
        }

        void _alignJoi(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size();
            NoGapAlignment::events_storage_t vec;
            vec.reserve(min(p_size, t_size) + 1);
            alignment_score_t score = 0, val; //, max_score = 0;

            vec.insert(vec.begin(), pattern[p_size - 1] != text[t_size - 1]);
            for (seq_len_t i = 1; i < min(p_size, t_size); ++i) {
                vec.insert(vec.begin(), pattern[p_size - i - 1] != text[t_size - i - 1]);
//                val = pattern[p_size - i - 1] == text[t_size - i - 1] ? _params.score.j_score.match : (_params.score.j_score.mism - _params.score.j_score.acc_mism*(pattern[p_size - i - 1] != text[t_size - i - 1]));
                score += pattern[p_size - i - 1] == text[t_size - i - 1] ? (_params.score.j_score.match + _params.score.j_score.acc_match * (pattern[p_size - i] == text[t_size - i])) : _params.score.j_score.mism;
                // max_score = std::max(max_score, score);
            }

            // if (max_score >= _params.threshold.j_threshold) {
            if (score >= _params.threshold.j_threshold) {
                avec->addAlignment(gene, p_size - min(t_size, p_size) + 1, t_size - min(t_size, p_size) + 1, vec);
            }
        }
        ///@}

    };


    //
    // Classic Smith-Waterman
    // Smith-Waterman aligner for finding maximal matches with gaps. Don't takes into account errors.
    //
    class SmithWatermanAligner : public VDJAlignerBase<GappedAlignmentVector> {
    public:

        SmithWatermanAligner()
        {
        }


        SmithWatermanAligner(const VDJRecombinationGenes &genes,
                             VDJAlignerParameters _params)
                : VDJAlignerBase<GappedAlignmentVector>(genes, _params)
        {
        }


    protected:

        /**
         *
         */
        ///@{
        void _alignVar(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, GappedAlignmentVector *avec) const {
            SWAlignmentMatrix mat(gene, pattern, text);

            mat.getBestAlignment(avec, pattern, text);
        }

        void _alignDiv(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, GappedAlignmentVector *avec) const {
            check_and_throw(false, "SWAlignerFunctor_D has not been implemented yet");
        }

        void _alignJoi(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, GappedAlignmentVector *avec) const {
            SWAlignmentMatrix mat(gene, pattern, text);

            mat.getBestAlignment(avec, pattern, text);
        }
        ///@}

    };


    //
    // Smith-Waterman with no gap allowed, but with errors
    // Smith-Waterman aligner without gaps, returns maximal matches with information about mismatch errors.
    //
    class SmithWatermanNoGapAligner : public VDJAlignerBase<ClonotypeNuc> {
    public:

        SmithWatermanNoGapAligner()
        {
        }


        SmithWatermanNoGapAligner(const VDJRecombinationGenes &genes,
                                  VDJAlignerParameters _params)
                : VDJAlignerBase<ClonotypeNuc>(genes, _params)
        {
        }


    protected:

        /**
         *
         */
        ///@{
        void _alignVar(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            SWNGAlignmentMatrix mat(gene, pattern, text);

            for (seq_len_t col_i = 0; col_i < text.size(); ++col_i) {
                for (seq_len_t row_i = 0; row_i < pattern.size(); ++row_i) {
                    mat.score(row_i + 1, col_i + 1) = std::max({mat.score(row_i, col_i) + (text[col_i] == pattern[row_i] ? _params.score.v_score.match : _params.score.v_score.mism), .0});
                }
            }

            mat.getBestAlignment(avec, pattern, text);
        }

        void _alignDiv(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            // avec->addAlignment(gene, );
            seq_len_t match_min_len = _params.min_D_len;
            seq_len_t t_size = text.size(), p_size = pattern.size(), min_size = min(t_size, p_size), min_subsize;
            seq_len_t p_start, t_start;
            AlignmentVectorBase::events_storage_t bitvec;
            bitvec.reserve(p_size + 1);

            for (seq_len_t pattern_i = 0; pattern_i < p_size - match_min_len + 1; ++pattern_i) {
                min_subsize = min(p_size - pattern_i, (int) t_size);
                if (min_subsize >= match_min_len) {
                    bitvec.resize(min_subsize);
                    for (seq_len_t i = 0; i < min_subsize; ++i) {
                        bitvec[i] = pattern[pattern_i + i] != text[i];
                    }
                    avec->addAlignment(gene, pattern_i + 1, 1, bitvec);
                }
            }

            for (seq_len_t text_i = 1; text_i < t_size - match_min_len + 1; ++text_i) {
                min_subsize = min((int) p_size, t_size - text_i);
                if (min_subsize >= match_min_len) {
                    bitvec.resize(min_subsize);
                    for (seq_len_t i = 0; i < min_subsize; ++i) {
                        bitvec[i] = pattern[i] != text[text_i + i];
                    }
                    avec->addAlignment(gene, 1, text_i + 1, bitvec);
                }
            }
        }

        void _alignJoi(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            SWNGAlignmentMatrix mat(gene, pattern, text);

            for (seq_len_t col_i = 0; col_i < text.size(); ++col_i) {
                for (seq_len_t row_i = 0; row_i < pattern.size(); ++row_i) {
                    mat.score(row_i + 1, col_i + 1) = std::max({mat.score(row_i, col_i) + (text[col_i] == pattern[row_i] ? _params.score.j_score.match : _params.score.j_score.mism), .0});
                }
            }

            mat.getBestAlignment(avec, pattern, text);
        }
        ///@}

    };
    

    // class NaiveAminoAcidAligner : public AbstractAligner {

    // public:

    //     NaiveAminoAcidAligner() { }


    //     virtual seq_len_t align5end(const string& pattern, const string& text) const {
    //         seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0, max_matches = 0, all_matches = 0;
    //         string codon_s = "";
    //         for (seq_len_t i = 0; i < std::min((seq_len_t) (1 + p_size / 3), t_size); ++i) {
    //             max_matches = 0;
    //             // go through all codons and find the maximal match
    //             CodonTable::Codons codon = _codons.codons(text[i]);
    //             while(!codon.end()) {
    //                 matches = 0;
    //                 codon_s = codon.next();
    //                 if (pattern[i*3] == codon_s[0] && i*3 < p_size) {
    //                     ++matches;
    //                     if (pattern[i*3 + 1] == codon_s[1] && (i*3 + 1) < p_size) {
    //                         ++matches;
    //                         if (pattern[i*3 + 2] == codon_s[2] && (i*3 + 2) < p_size) {
    //                             ++matches;
    //                         }
    //                     }
    //                 }

    //                 if (matches > max_matches) {
    //                     max_matches = matches;
    //                     if (max_matches == 3) { break; }
    //                 }
    //             }

    //             // if match == 3 then go to the next amino acid
    //             all_matches += max_matches;
    //             if (max_matches != 3) { break; }
    //         }

    //         return all_matches;
    //     }


    //     virtual seq_len_t align3end(const string& pattern, const string& text) const {
    //         seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0, max_matches = 0, all_matches = 0;
    //         string codon_s = "";
    //         for (seq_len_t i = 0; i < std::min((seq_len_t) (1 + p_size / 3), t_size); ++i) {
    //             max_matches = 0;
    //             // go through all codons and find the maximal match
    //             CodonTable::Codons codon = _codons.codons(text[t_size - i - 1]);
    //             while(!codon.end()) {
    //                 matches = 0;
    //                 codon_s = codon.next();
    //                 if (pattern[p_size - 1 - i*3] == codon_s[2] && i*3 < p_size) {
    //                     ++matches;
    //                     if (pattern[p_size - 1 - i*3 - 1] == codon_s[1] && (i*3 + 1) < p_size) {
    //                         ++matches;
    //                         if (pattern[p_size - 1 - i*3 - 2] == codon_s[0] && (i*3 + 2) < p_size) {
    //                             ++matches;
    //                         }
    //                     }
    //                 }

    //                 if (matches > max_matches) {
    //                     max_matches = matches;
    //                     if (max_matches == 3) { break; }
    //                 }
    //             }

    //             // if match == 3 then go to the next amino acid
    //             all_matches += max_matches;
    //             if (max_matches != 3) { break; }
    //         }

    //         return all_matches;
    //     }


    //     virtual LocalAlignmentIndices alignLocal(const string& pattern, const string& text, seq_len_t match_min_len = 3) const {
    //         return LocalAlignmentIndices(vector<seq_len_t>(1));
    //     }

    // protected:

    //     const CodonTable _codons;

    // };


}

#endif