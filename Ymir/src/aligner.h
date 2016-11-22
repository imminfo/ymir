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
#include "repertoire.h"


using namespace std;


namespace ymir {

    template <typename ClonotypeType>
    class VDJAlignerBase;

    class NaiveCDR3NucleotideAligner;

    class CDR3NucleotideAligner;

    class SmithWatermanAligner;

    class SmithWatermanNoGapAligner;

    class CDR3AminoAcidAligner;


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
            this->_alignVar(id, _genes.V()[id].sequence, sequence, &vec);
            return vec;
        }

        alignment_vector_t alignDiv(seg_index_t id, const sequence_t &sequence) const {
            alignment_vector_t vec;
            this->_alignDiv(id, _genes.D()[id].sequence, sequence, &vec);
            return vec;
        }

        alignment_vector_t alignJoi(seg_index_t id, const sequence_t &sequence) const {
            alignment_vector_t vec;
            this->_alignJoi(id, _genes.J()[id].sequence, sequence, &vec);
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
                this->_alignVar(id, _genes.V()[id].sequence, this->sequence(), &vec);
            }
            this->addVarAlignment(vec);
            return vec.size() != 0;
        }

        bool alignVar(seg_index_t id) {
            alignment_vector_t vec;
            this->_alignVar(id, _genes.V()[id].sequence, this->sequence(), &vec);
            this->addVarAlignment(vec);
            return vec.size() != 0;
        }

        bool alignDiv() {
            alignment_vector_t vec;
            for (seg_index_t id = 1; id <= _genes.D().max(); ++id) {
                vec.clear(); // TODO: some strange behaviour here after I added this line. Be careful. Maybe there is some bug in clear().
                this->_alignDiv(id, _genes.D()[id].sequence, this->sequence(), &vec);
                if (vec.size()) {
                    this->addDivAlignment(vec);
                }
            }
            return vec.size() != 0;
        }

        bool alignDiv(seg_index_t id) {
            alignment_vector_t vec;
            this->_alignDiv(id, _genes.D()[id].sequence, this->sequence(), &vec);
            this->addDivAlignment(vec);
            return vec.size() != 0;
        }

        bool alignJoi() {
            alignment_vector_t vec;
            for (seg_index_t id = 1; id <= _genes.J().max(); ++id) {
                this->_alignJoi(id, _genes.J()[id].sequence, this->sequence(), &vec);
            }
            this->addJoiAlignment(vec);
            return vec.size() != 0;
        }

        bool alignJoi(seg_index_t id) {
            alignment_vector_t vec;
            this->_alignJoi(id, _genes.J()[id].sequence, this->sequence(), &vec);
            this->addJoiAlignment(vec);
            return vec.size() != 0;
        }

        bool align(GeneSegments gene_segment, seg_index_t id) {
            switch (gene_segment) {
                case VARIABLE:
                    return this->alignVar(id);
                case JOINING:
                    return this->alignJoi(id);
                case DIVERSITY:
                    return this->alignDiv(id);
                default:
                    return false;
            }
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


    typedef VDJAlignerBase<ClonotypeNuc> VDJNucAligner;


    typedef VDJAlignerBase<ClonotypeAA> VDJAAAligner;


    /**
     * \class NaiveCDR3NucleotideAligner
     *
     * \brief CDR3-only alignment without errors - align V starting from the left edge,
     * J starting from the right edge, and align D everywhere.
     */
    class NaiveCDR3NucleotideAligner : public VDJNucAligner {
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
    class CDR3NucleotideAligner : public VDJNucAligner {
    public:

        CDR3NucleotideAligner()
        {
        }


        CDR3NucleotideAligner(const VDJRecombinationGenes &genes,
                              VDJAlignerParameters params)
                : VDJNucAligner(genes, params)
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
//    class SmithWatermanAligner : public VDJAlignerBase<GappedAlignmentVector> {
    class SmithWatermanAligner : public VDJNucAligner {
    public:

        SmithWatermanAligner()
        {
        }


        SmithWatermanAligner(const VDJRecombinationGenes &genes,
                             VDJAlignerParameters params)
//                : VDJAlignerBase<GappedAlignmentVector>(genes, _params)
        : VDJNucAligner(genes, _params)
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
    class SmithWatermanNoGapAligner : public VDJNucAligner {
    public:

        SmithWatermanNoGapAligner()
        {
        }


        SmithWatermanNoGapAligner(const VDJRecombinationGenes &genes,
                                  VDJAlignerParameters params)
                : VDJNucAligner(genes, params)
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


    //
    // CDR3-only alignment for amino acid sequences - align V starting from the left edge,
    // J starting from the right edge, and align D everywhere.
    //
    class CDR3AminoAcidAligner : public VDJAAAligner {
    public:

        CDR3AminoAcidAligner()
        {
        }


        CDR3AminoAcidAligner(const VDJRecombinationGenes &genes,
                              VDJAlignerParameters params)
                : VDJAAAligner(genes, params)
        {
        }


        ClonotypeAA toAminoAcid(const ClonotypeNuc &clonotype, bool use_aligned_segments = true) {
            this->setSequence(clonotype.aa_sequence());
            this->setRecombination(clonotype.recombination());
            if (use_aligned_segments) {
                for (seg_index_t id = 0; id < clonotype.nVar(); ++id) {
                    this->alignVar(clonotype.getVar(id));
                }
                for (seg_index_t id = 0; id < clonotype.nJoi(); ++id) {
                    this->alignJoi(clonotype.getJoi(id));
                }
            } else {
                this->alignVar();
                this->alignJoi();
            }
            if (_recomb == VDJ_RECOMB) { this->alignDiv(); }

            return this->buildClonotype();
        }


        std::vector<size_t> toAminoAcid(const ClonesetViewNuc &cloneset, ClonesetAA *cloneset_aa, bool use_aligned_segments = true) {
            ClonotypeAAVector vec;
            vec.reserve(cloneset.size());
            std::vector<size_t> noncoding_indices;
            noncoding_indices.reserve(cloneset.size());
            size_t stats_bad_v = 0, stats_bad_j = 0, stats_bad_both = 0;
            for (size_t i = 0; i < cloneset.size(); ++i) {
                if (cloneset[i].isCoding()) {
                    auto res = this->toAminoAcid(cloneset[i], use_aligned_segments);
                    if (res.nVar() && res.nJoi()) {
                        vec.push_back(res);
                    } else {
                        stats_bad_v    += (!res.nVar()) &&   res.nJoi();
                        stats_bad_j    +=   res.nVar()  && (!res.nJoi());
                        stats_bad_both += (!res.nVar()) && (!res.nJoi());
                    }
                } else {
                    noncoding_indices.push_back(i);
                }

                if ((i+1) % 50000 == 0) {
                    std::cout << "converted " << (size_t) (i+1) << " clonotypes" << std::endl;
                }
            }
            std::cout << "err. clonotypes: " << std::endl <<
                    "\toverall: " << (size_t) (stats_bad_v + stats_bad_j + stats_bad_both) << std::endl <<
                    "\tbad V:   "   << (size_t) (stats_bad_v) << std::endl <<
                    "\tbad J:   "   << (size_t) (stats_bad_j) << std::endl <<
                    "\tbad V+J: " << (size_t) (stats_bad_both) << std::endl;
            cloneset_aa->swap(vec);

            return noncoding_indices;
        }


    protected:

        /**
         *
         * For Diversity gene each inner codon's position corresponds to possible codons
         * on the corresponding position on the D gene sequence, i.e., without any modification
         * from neighbor positions. But each outer positions undergone modifications same as for
         * V / J sequences. For the D3' side i-th position depends on i+1-th position.
         * For the D5' side i+1-th position depends on i-th position.
         *
         */
        ///@{
        void _alignVar(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, CodonAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0, max_matches = 0, all_matches = 0;
            CodonAlignmentVector::events_storage_t bits;
            bool is_ok, stop = false;
            size_t bits_size;
            int num_codons;

            int max_iter = std::min((seq_len_t) (p_size / 3 + static_cast<seq_len_t>((p_size % 3) != 0)), t_size);
            for (seq_len_t i = 0; i < max_iter; ++i) {
                num_codons = 0;
                is_ok = false;

                // go through all codons and find the maximal match
                if (i*3 < p_size) {
                    is_ok = CodonTable::table().check_nucl(text[i], pattern[i*3], 0, &bits);
                    ++num_codons;
                }

                if (is_ok && (i*3 + 1) < p_size) {
                    is_ok = CodonTable::table().check_nucl(text[i], pattern[i*3 + 1], 1, &bits);
                    ++num_codons;
                } else {
                    is_ok = false;
                }

                if (is_ok) {
                    // find intersected codons
                    for (int bit_i = 1; bit_i <= 6; ++bit_i) {
                        bits[bits.size() - bit_i] = bits[bits.size() - bit_i] & bits[bits.size() - 6 - bit_i];
                    }

                    if ((i*3 + 2) < p_size) {
                        is_ok = CodonTable::table().check_nucl(text[i], pattern[i*3 + 2], 2, &bits);
                        ++num_codons;

                        if (is_ok) {
                            // find intersected codons
                            for (int bit_i = 1; bit_i <= 6; ++bit_i) {
                                bits[bits.size() - bit_i] = bits[bits.size() - bit_i] & bits[bits.size() - 6 - bit_i];
                            }
                        }
                    }
                }

//                std::cout << (int) num_codons << std::endl;

                // remove trailing zeros
                for (int codon_i = 0; codon_i < num_codons; ++codon_i) {
                    int bitsum = 0;
                    for (int bit_i = 1; bit_i <= 6; ++bit_i) {
                        bitsum += bits[bits.size() - bit_i];
                    }
                    if (!bitsum) {
                        bits.resize(bits.size() - 6);
                        stop = true;
                    } else {
                        break;
                    }
                }

                if (stop) { break; }
            }

            if (bits.size() / 6 >= _params.threshold.v_threshold) { avec->addAlignment(gene, 1, 1, bits); }
//            avec->addAlignment(gene, 1, 1, bits);
        }

        void _alignDiv(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, CodonAlignmentVector *avec) const {
            seq_len_t match_min_len = _params.min_D_len;
            seq_len_t t_size = text.size() * 3, p_size = pattern.size(), min_size = min(t_size, p_size), min_subsize;
            bool open_match;
            seq_len_t p_start, t_start;
            CodonAlignmentVector bits;
            std::vector<bool> nucbits;
            bool is_ok;
            int bitsum;

            for (seq_len_t pattern_i = 0; pattern_i < p_size - match_min_len + 1; ++pattern_i) {
                nucbits.clear();

                open_match = false;
                min_subsize = min(p_size - pattern_i, (int) t_size);
                for (seq_len_t i = 0; i < min_subsize; ++i) {
                    is_ok = CodonTable::table().check_nucl(text[i / 3], pattern[pattern_i + i], i % 3, &nucbits);

                    // find intersected codons, if the i-th position isn't the first one in the processed codon
                    if ((nucbits.size() > 6) && ((i % 3) > 0)) {
                        bitsum = 0;
                        for (int bit_i = 1; bit_i <= 6; ++bit_i) {
                            bitsum += nucbits[nucbits.size() - bit_i] & nucbits[nucbits.size() - 6 - bit_i];
                        }
                        is_ok &= static_cast<bool>(bitsum);
                    }

                    if (is_ok) {
                        if (!open_match) {
                            p_start = pattern_i + i;
                            t_start = i;
                            open_match = true;
                        }
                    } else {
                        if (open_match) {
                            if (nucbits.size()) { nucbits.resize(nucbits.size() - 6); }
                            if ((pattern_i + i - p_start) >= match_min_len) {
                                // Save each of the 6 codon bits for each nucleotide position
                                avec->addAlignment(gene, p_start + 1, t_start + 1, nucbits);
                            }
                            open_match = false;
                        }
                        nucbits.clear();
                    }
                }
                if (open_match && (pattern_i + min_subsize - p_start) >= match_min_len) {
                    avec->addAlignment(gene, p_start + 1, t_start + 1, nucbits);
                    nucbits.clear();
                }
            }

            for (seq_len_t text_i = 1; text_i < t_size - match_min_len + 1; ++text_i) {
                nucbits.clear();

                open_match = false;
                min_subsize = min((int) p_size, t_size - text_i);
                for (seq_len_t i = 0; i < min_subsize; ++i) {
                    is_ok = CodonTable::table().check_nucl(text[(text_i + i) / 3], pattern[i], (text_i + i) % 3, &nucbits);

                    if ((nucbits.size() > 6) && (((text_i + i) % 3) > 0)) {
                        bitsum = 0;
                        for (int bit_i = 1; bit_i <= 6; ++bit_i) {
                            bitsum += nucbits[nucbits.size() - bit_i] & nucbits[nucbits.size() - 6 - bit_i];
                        }
                        is_ok &= static_cast<bool>(bitsum);
                    }

                    if (is_ok) {
                        if (!open_match) {
                            p_start = i;
                            t_start = text_i + i;
                            open_match = true;
                        }
                    } else {
                        if (open_match) {
                            if (nucbits.size()) { nucbits.resize(nucbits.size() - 6); }
                            if ((i - p_start) >= match_min_len) {
                                avec->addAlignment(gene, p_start + 1, t_start + 1, nucbits);
                            }
                            open_match = false;
                        }
                        nucbits.clear();
                    }
                }
                if (open_match && (min_subsize - p_start) >= match_min_len) {
                    avec->addAlignment(gene, p_start + 1, t_start + 1, nucbits);
                    nucbits.clear();
                }
            }
        }

        void _alignJoi(seg_index_t gene, const sequence_t &pattern, const sequence_t &text, CodonAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size();
            CodonAlignmentVector::events_storage_t bits;
            bool is_ok, stop = false;
            size_t bits_size;
            int num_codons, all_codons = 0;

            int max_iter = std::min((seq_len_t) (p_size / 3 + static_cast<seq_len_t>((p_size % 3) != 0)), t_size);
            for (seq_len_t i = 0; i < max_iter; ++i) {
                num_codons = 0;
                is_ok = false;

                // go through all codons and find the maximal match
                if (i*3 < p_size) {
                    is_ok = CodonTable::table().check_nucl(text[t_size - 1 - i], pattern[p_size - 1 - i*3], 2, &bits);
                    ++num_codons;
                }

                if (is_ok && (i*3 + 1) < p_size) {
                    is_ok = CodonTable::table().check_nucl(text[t_size - 1 - i], pattern[p_size - 2 - i*3], 1, &bits);
                    ++num_codons;
                } else {
                    is_ok = false;
                }

                if (is_ok) {
                    // find intersected codons
                    for (int bit_i = 1; bit_i <= 6; ++bit_i) {
                        bits[bits.size() - bit_i] = bits[bits.size() - bit_i] & bits[bits.size() - 6 - bit_i];
                    }

                    if ((i*3 + 2) < p_size) {
                        is_ok = CodonTable::table().check_nucl(text[t_size - 1 - i], pattern[p_size - 3 - i*3], 0, &bits);
                        ++num_codons;

                        if (is_ok) {
                            // find intersected codons
                            for (int bit_i = 1; bit_i <= 6; ++bit_i) {
                                bits[bits.size() - bit_i] = bits[bits.size() - bit_i] & bits[bits.size() - 6 - bit_i];
                            }
                        }
                    }
                }

                // remove trailing zeros
                all_codons += num_codons;
                for (int codon_i = 0; codon_i < num_codons; ++codon_i) {
                    int bitsum = 0;
                    for (int bit_i = 1; bit_i <= 6; ++bit_i) {
                        bitsum += bits[bits.size() - bit_i];
                    }
                    if (!bitsum) {
                        bits.resize(bits.size() - 6);
                        --all_codons;
                        stop = true;
                    } else {
                        break;
                    }
                }

                if (stop) { break; }
            }

            // reverse bits
            CodonAlignmentVector::events_storage_t bits2;
            bits2.resize(bits.size());
            for (size_t bit_i = 0; bit_i < bits.size(); bit_i += 6) {
                bits2[bit_i] = bits[bits.size() - bit_i - 6];
                bits2[bit_i + 1] = bits[bits.size() - bit_i - 5];
                bits2[bit_i + 2] = bits[bits.size() - bit_i - 4];
                bits2[bit_i + 3] = bits[bits.size() - bit_i - 3];
                bits2[bit_i + 4] = bits[bits.size() - bit_i - 2];
                bits2[bit_i + 5] = bits[bits.size() - bit_i - 1];
            }

            if (bits2.size() / 6 >= _params.threshold.j_threshold) { avec->addAlignment(gene, p_size + 1 - bits2.size() / 6, t_size*3 + 1 - all_codons, bits2); }
//            avec->addAlignment(gene, p_size + 1 - bits2.size() / 6, t_size*3 + 1 - all_codons, bits2);
        }
        ///@}

    };
    

    // class NaiveAminoAcidAligner : public AbstractAligner {

    // public:

    //     NaiveAminoAcidAligner() { }


    //     virtual seq_len_t align5end(const string& pattern, const string& text) const {
//             seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0, max_matches = 0, all_matches = 0;
//             string codon_s = "";
//             for (seq_len_t i = 0; i < std::min((seq_len_t) (1 + p_size / 3), t_size); ++i) {
//                 max_matches = 0;
//                 // go through all codons and find the maximal match
//                 CodonTable::Codons codon = _codons.codons(text[i]);
//                 while(!codon.end()) {
//                     matches = 0;
//                     codon_s = codon.next();
//                     if (pattern[i*3] == codon_s[0] && i*3 < p_size) {
//                         ++matches;
//                         if (pattern[i*3 + 1] == codon_s[1] && (i*3 + 1) < p_size) {
//                             ++matches;
//                             if (pattern[i*3 + 2] == codon_s[2] && (i*3 + 2) < p_size) {
//                                 ++matches;
//                             }
//                         }
//                     }
//
//                     if (matches > max_matches) {
//                         max_matches = matches;
//                         if (max_matches == 3) { break; }
//                     }
//                 }
//
//                 // if match == 3 then go to the next amino acid
//                 all_matches += max_matches;
//                 if (max_matches != 3) { break; }
//             }
//
//             return all_matches;
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