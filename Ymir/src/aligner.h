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


#include "alignment.h"
#include "genesegment.h"


using namespace std;


namespace ymir {


    /**
     * \class VDJAlignerBase
     */
    // template <typename A, typename V_Aligner, typename D_Aligner, typename J_Aligner>
    // class VDJAlignerBase {
    // public:

    //     /**
    //      * \struct VDJAlignment
    //      */
    //     struct VDJAlignment {


    //         VDJAlignment()
    //                 : _gene(UNDEF_GENE), _n_alignments(0), _alignments(nullptr)
    //         { }


    //         VDJAlignment(GeneSegments gene, seg_index_t n_alignments, A *alignments)
    //                 : _gene(gene), _n_alignments(n_alignments), _alignments(alignments)
    //         { }


    //         VDJAlignment(const VDJAlignment &other)
    //                 : _gene(other._gene), _n_alignments(other._n_alignments)
    //         {
    //             _alignments = new SegmentAlignment[other._n_alignments];
    //             std::copy(other._alignments, other._alignments + other._n_alignments, _alignments);
    //         }


    //         ~VDJAlignment() {
    //             delete [] _alignments;
    //         }


    //         VDJAlignment& operator=(const VDJAlignment &other) {
    //             _gene = other._gene;
    //             _n_alignments = other._n_alignments;
    //             _alignments = new SegmentAlignment[other._n_alignments];
    //             std::copy(other._alignments, other._alignments + other._n_alignments, _alignments);
    //         }


    //         GeneSegments gene() const { return _gene; }


    //         size_t size() const { return _n_alignments; }


    //         const A& getAlignment(seg_index_t i) const { return _alignments[i]; }

    //     protected:

    //         GeneSegments _gene;
    //         A *_alignments;
    //         seg_index_t _n_alignments;

    //     };


    //     VDJAlignerBase() { }


    //     /**
    //      *
    //      */
    //     VDJAlignerBase(const VDJRecombinationGenes &genes,
    //                        alignment_score_t threshold
    //                        // const AlignmentEventScore &v_score = AlignmentEventScore(1, -1, -1, -1),
    //                        // const AlignmentEventScore &d_score = AlignmentEventScore(1, -1, -1, -1),
    //                        // const AlignmentEventScore &j_score = AlignmentEventScore(1, -1, -1, -1))
    //             : _genes(genes), _threshold(threshold)
    //     {
    //     }


    //     ~VDJAlignerBase() {
    //     }


    //     /**
    //      * \brief Align the given sequence to the specific gene segment of the specific gene.
    //      *
    //      * \param sequence Pattern sequence.
    //      * \param seg_index Index of the target gene segment.
    //      */
    //     ///@{
    //     SegmentAlignment alignVar(const sequence_t &sequence, seg_index_t seg_index, const sequence_t &segment_seq) const {
    //         this->alignOneSegment<V_Aligner>(sequence, seg_index, segment_seq);
    //     }

    //     SegmentAlignmentVector alignDiv(const sequence_t &sequence, seg_index_t seg_index, const sequence_t &segment_seq) const {
    //         this->alignManySegments<D_Aligner>(sequence, seg_index, segment_seq);
    //     }

    //     SegmentAlignment alignJoi(const sequence_t &sequence, seg_index_t seg_index, const sequence_t &segment_seq) const {
    //         this->alignOneSegment<J_Aligner>(sequence, seg_index, segment_seq);
    //     }
    //     ///@}


    //     /**
    //      * \brief Align the given sequence to all gene segments of the specific gene.
    //      *
    //      * \param sequence Pattern sequence.
    //      */
    //     ///@{
    //     void alignVar(const sequence_t& sequence) {
    //         this->alignGeneSegments<V_Aligner>(sequence, _genes.V(), _last_v_alignment);
    //     }

    //     void alignDiv(const sequence_t &sequence) {
    //         this->alignGeneSegments<D_Aligner>(sequence, _genes.D(), _last_d_alignment);
    //     }

    //     void alignJoi(const sequence_t &sequence) {
    //         this->alignGeneSegments<J_Aligner>(sequence, _genes.J(), _last_j_alignment);
    //     }
    //     ///@}


    //     /**
    //      * \brief Access the latest alignment results.
    //      */
    //     ///@{
    //     const VDJAlignment& lastVarAlignment() const { return _last_v_alignment; }

    //     const VDJAlignment& lastDivAlignment() const { return _last_d_alignment; }

    //     const VDJAlignment& lastJoiAlignment() const { return _last_j_alignment; }
    //     ///@}

    // protected:

    //     alignment_score_t _threshold;
    //     VDJRecombinationGenes _genes;
    //     VDJAlignment _last_v_alignment, _last_d_alignment, _last_j_alignment;


    //     template <typename F_Aligner>
    //     SegmentAlignment alignOneSegment(const sequence_t &sequence, seg_index_t seg_index, const sequence_t &segment_seq) const {
    //         return F_Aligner(sequence, seg_index, segment_seq);
    //     }


    //     template <typename F_Aligner>
    //     void alignGeneSegments(const sequence_t& sequence, const GeneSegmentAlphabet &gsa, VDJAlignment &gs_alignment) {
    //         SegmentAlignmentVector vec;
    //         vec.reserve(gsa.size() / 2 + 2);
    //         SegmentAlignment tmp;
    //         for (seg_index_t i = 1; i <= gsa.max(); ++i) {
    //             tmp = F_Aligner(sequence, i, gsa[i].sequence);
    //             if (tmp.score >= _threshold) {
    //                 vec.push_back(tmp);
    //             }
    //         }
    //         SegmentAlignment *arr = new SegmentAlignment[vec.size()];
    //         gs_alignment = VDJAlignment(gsa.gene_segment(), vec.size(), arr);
    //     }

    //     template <typename F_Aligner>
    //     void alignManySegments(const sequence_t& sequence, const GeneSegmentAlphabet &gsa, VDJAlignment &gs_alignment) {
    //         SegmentAlignmentVector vec, tmp;
    //         vec.reserve(gsa.size() * 4);
    //         for (seg_index_t i = 1; i <= gsa.max(); ++i) {
    //             tmp = F_Aligner(sequence, i, gsa[i].sequence);
    //             vec.insert(vec.end(), tmp.begin(), tmp.end());
    //         }
    //         SegmentAlignment *arr = new SegmentAlignment[vec.size()];
    //         gs_alignment = VDJAlignment(gsa.gene_segment(), vec.size(), arr);
    //     }

    // };


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


    //
    // CDR3-only alignment without errors - align V starting from the left edge, 
    // J starting from the right edge, and align D everywhere.
    //

    /**
     * \brief
     */
    ///@{
    struct NaiveCDR3AlignerFunctor_V {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0;

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[i] != text[i]) { break; }
                matches += 1;
            }

            avec->addAlignment(0, 1, 1, matches);
        }
    };

    struct NaiveCDR3AlignerFunctor_D {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec, seq_len_t match_min_len) const {
            seq_len_t t_size = text.size(), p_size = pattern.size(), min_size = min(t_size, p_size), min_subsize;
            bool open_match;
            seq_len_t p_start, t_start;

            for (seq_len_t pattern_i = 1 /* WTF?! */ ; pattern_i < p_size - match_min_len + 1; ++pattern_i) {
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
                            avec->addAlignment(0, p_start + 1, t_start + 1, pattern_i + i - p_start);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (pattern_i + min_subsize - p_start) >= match_min_len) {
                    avec->addAlignment(0, p_start + 1, t_start + 1, pattern_i + min_subsize - p_start);
                }
            }

            for (seq_len_t text_i = 0; text_i < t_size - match_min_len + 1; ++text_i) {
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
                            avec->addAlignment(0, p_start + 1, t_start + 1, i - p_start);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (min_subsize - p_start) >= match_min_len) {
                    avec->addAlignment(0, p_start + 1, t_start + 1, min_subsize - p_start);
                }
            }
        }
    };

    struct NaiveCDR3AlignerFunctor_J {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0;

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[p_size - i - 1] != text[t_size - i - 1]) { break; }
                matches += 1;
            }

            avec->addAlignment(0, 1, t_size - matches + 1, matches);
        }
    };
    ///@}


    // typedef VDJAlignerBase<NaiveCDR3AlignerFunctor_V, NaiveCDR3AlignerFunctor_D, NaiveCDR3AlignerFunctor_J> NaiveCDR3NucleotideAligner;


    //
    // CDR3-only alignment with errors - align V starting from the left edge, 
    // J starting from the right edge, and align D everywhere.
    //

    /**
     * \brief
     */
    ///@{
    struct CDR3AlignerFunctor_V {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size();
            NoGapAlignment::events_storage_t vec;
            vec.reserve(min(p_size, t_size));

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[i] == text[i]) { 
                    vec.push_back(false);
                } else {
                    vec.push_back(true);
                }
            }

            avec->addAlignment(0, 1, 1, vec);
        }
    };

    struct CDR3AlignerFunctor_D {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec, seq_len_t match_min_len) const {
            // pass
        }
    };

    struct CDR3AlignerFunctor_J {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            seq_len_t p_size = pattern.size(), t_size = text.size();
            NoGapAlignment::events_storage_t vec;
            vec.reserve(min(p_size, t_size));

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[p_size - i - 1] != text[t_size - i - 1]) {
                    vec.push_back(false);
                } else {
                    vec.push_back(true);
                }
            }

            avec->addAlignment(0, 0, 1, t_size - min(t_size, p_size) + 1, vec);
        }
    };
    ///@}


    // typedef VDJAlignerBase<CDR3AlignerFunctor_V, CDR3AlignerFunctor_D, CDR3AlignerFunctor_J> CDR3NucleotideAligner;


    // /**
    //  * \class AlignmentMatrix
    //  */
    // class AlignmentMatrix {
    // public:


    //     AlignmentMatrix(seq_len_t nrow, seq_len_t ncol)
    //             : _nrow(nrow), _ncol(ncol), _events((nrow+1) * (ncol + 1))
    //     {
    //         _starts = new bool[nrow * ncol];
    //         std::fill(_starts, _starts + nrow * ncol, false);
    //     }


    //     ~AlignmentMatrix() {
    //         delete [] _starts;
    //     }


    //     // alignment_event_t getEvent(seq_len_t row, seq_len_t col) const { return _events[row * _nrow + col]; }


    //     void setEvent(seq_len_t row, seq_len_t col, const alignment_event_t &event) { _events.setEvent(row * _nrow + col, event); }


    //     bool is_start(seq_len_t row, seq_len_t col) const { return _starts[row * _nrow + col]; }


    // private:

    //     seq_len_t _nrow, _ncol;
    //     AlignmentEventVector _events;
    //     bool *_starts;

    // };


    //
    // Classic Smith-Waterman
    //

    /**
     * \brief
     */
    ///@{
    struct SWAlignerFunctor_VJ {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            // avec->addAlignment(0, );
        }
    };

    struct SWAlignerFunctor_D {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            // avec->addAlignment(0, );
        }
    };
    ///@}


    /**
     * \typedef SmithWatermanAligner
     *
     * \brief Smith-Waterman aligner for finding maximal matches with gaps. Don't takes into account errors.
     */
    // typedef VDJAlignerBase<SWAlignerFunctor_VJ, SWAlignerFunctor_D, SWAlignerFunctor_VJ> SmithWatermanAligner;


    //
    // Smith-Waterman with no gap allowed, but with errors
    //

    /**
     * \brief
     */
    ///@{
    struct SWNGAlignerFunctor_VJ {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            // avec->addAlignment(0, );
        }
    };

    struct SWNGAlignerFunctor_D {
        void operator()(const sequence_t &pattern, const sequence_t &text, NoGapAlignmentVector *avec) const {
            // avec->addAlignment(0, );
        }
    };
    ///@}


    /**
     * \typedef SmithWatermanNoGapAligner
     *
     * \brief Smith-Waterman aligner without gaps, returns maximal matches with information about mismatch errors.
     */
    // typedef VDJAlignerBase<SWNGAlignerFunctor_VJ, SWNGAlignerFunctor_D, SWNGAlignerFunctor_VJ> SmithWatermanNoGapAligner;

}

#endif