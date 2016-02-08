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


#include "vdj_alignment_builder.h"
#include "genesegment.h"


using namespace std;


namespace ymir {

    /**
     * \struct VDJAlignerParameters
     */
    struct VDJAlignerParameters {


        static const alignment_score_t default_thr = 1;


        static const seq_len_t default_minlen = 3;


        VDJAlignerParameters() 
            : threshold(default_thr), 
              min_D_len(default_minlen)
        {
        }


        VDJAlignerParameters(alignment_score_t thr, seq_len_t minlen) 
            : threshold(thr), 
              min_D_len(minlen)
        {                
        }


        alignment_score_t threshold; ///
        seq_len_t min_D_len; ///

    private:

    };


    /**
     * \class VDJAlignerBase
     */
    template <typename AlignmentType, typename V_Aligner, typename D_Aligner, typename J_Aligner>
    class VDJAlignerBase {

    public:

        /**
         *
         */
        VDJAlignerBase(const VDJRecombinationGenes &genes,
                       VDJAlignerParameters params)
                           // const AlignmentEventScore &v_score = AlignmentEventScore(1, -1, -1, -1),
                           // const AlignmentEventScore &d_score = AlignmentEventScore(1, -1, -1, -1),
                           // const AlignmentEventScore &j_score = AlignmentEventScore(1, -1, -1, -1))
                : _genes(genes), _params(params)
        {
        }


        ~VDJAlignerBase() 
        {
        }


        /**
         * \brief General methods for alignment.
         */
        ///@{
        AlignmentType alignVar(seg_index_t id, const sequence_t &pattern) const {
            AlignmentType vec;
            this->_alignVar(id, pattern, &vec);
            return vec;
        }

        AlignmentType alignDiv(seg_index_t id, const sequence_t &pattern) const {
            AlignmentType vec;
            this->_alignDiv(id, pattern, &vec);
            return vec;
        }

        AlignmentType alignJoi(seg_index_t id, const sequence_t &pattern) const {
            AlignmentType vec;
            this->_alignJoi(id, pattern, &vec);
            return vec;
        }
        ///@}


        /**
         * \brief Align the given sequence to all gene segments of the specific gene.
         *
         * \param sequence Pattern sequence.
         */
        ///@{
        void alignVar(const sequence_t &sequence) {
            NoGapAlignmentVector vec;
            for (seg_index_t id = 1; id <= _genes.V().max(); ++id) {
                this->_alignVar(id, sequence, &vec);
            }
            _builder.addVarAlignment(vec);
        }

        void alignDiv(const sequence_t &sequence) {
            NoGapAlignmentVector vec;
            for (seg_index_t id = 1; id <= _genes.D().max(); ++id) {
                this->_alignDiv(id, sequence, &vec);
                _builder.addDivAlignment(vec);
            }
        }

        void alignJoi(const sequence_t &sequence) {
            NoGapAlignmentVector vec;
            for (seg_index_t id = 1; id <= _genes.J().max(); ++id) {
                this->_alignJoi(id, sequence, &vec);
            }
            _builder.addJoiAlignment(vec);
        }
        ///@}


        /**
         * \brief Access the latest alignment results and optionally move it.
         */
         VDJAlignment&& last_alignment_and_clear() { return _builder.build(); }

    protected:

        VDJAlignerParameters _params;
        VDJRecombinationGenes _genes;
        VDJAlignmentBuilder _builder;

        V_Aligner _V_Aligner;
        D_Aligner _D_Aligner;
        J_Aligner _J_Aligner;


        VDJAlignerBase()
        {
        }


        /**
         * \brief Internal alignment functions.
         *
         * \param id
         * \param pattern
         * \param vec
         */
        ///@{
        void _alignVar(seg_index_t id, const sequence_t &pattern, AlignmentType *vec) const {
            _V_Aligner(id, pattern, _genes.V()[id].sequence, vec, _params);
        }

        void _alignDiv(seg_index_t id, const sequence_t &pattern, AlignmentType *vec) const {
            _D_Aligner(id, pattern, _genes.D()[id].sequence, vec, _params);
        }

        void _alignJoi(seg_index_t id, const sequence_t &pattern, AlignmentType *vec) const {
            _J_Aligner(id, pattern, _genes.J()[id].sequence, vec, _params);
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


    //
    // CDR3-only alignment without errors - align V starting from the left edge, 
    // J starting from the right edge, and align D everywhere.
    //


    /**
     * \brief
     */
    ///@{
    struct NaiveCDR3AlignerFunctor_V {
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec, 
                        VDJAlignerParameters params = VDJAlignerParameters()) const 
        {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0;

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[i] != text[i]) { break; }
                matches += 1;
            }

            avec->addAlignment(gene, 1, 1, matches);
        }
    };

    struct NaiveCDR3AlignerFunctor_D {
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec, 
                        VDJAlignerParameters params = VDJAlignerParameters()) const 
        {
            seq_len_t match_min_len = params.min_D_len;
            seq_len_t t_size = text.size(), p_size = pattern.size(), min_size = min(t_size, p_size), min_subsize;
            bool open_match;
            seq_len_t p_start, t_start;

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
                            avec->addAlignment(gene, p_start + 1, t_start + 1, i - p_start);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (min_subsize - p_start) >= match_min_len) {
                    avec->addAlignment(gene, p_start + 1, t_start + 1, min_subsize - p_start);
                }
            }

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
                            avec->addAlignment(gene, p_start + 1, t_start + 1, pattern_i + i - p_start);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (pattern_i + min_subsize - p_start) >= match_min_len) {
                    avec->addAlignment(gene, p_start + 1, t_start + 1, pattern_i + min_subsize - p_start);
                }
            }
        }
    };

    struct NaiveCDR3AlignerFunctor_J {
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec,
                        VDJAlignerParameters params = VDJAlignerParameters()) const
        {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0;

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[p_size - i - 1] != text[t_size - i - 1]) { break; }
                matches += 1;
            }

            avec->addAlignment(gene, p_size - matches + 1, t_size - matches + 1, matches);
        }
    };
    ///@}


    /**
     *
     */
    typedef VDJAlignerBase<NoGapAlignmentVector,
                           NaiveCDR3AlignerFunctor_V, 
                           NaiveCDR3AlignerFunctor_D, 
                           NaiveCDR3AlignerFunctor_J> NaiveCDR3NucleotideAligner;


    //
    // CDR3-only alignment with errors - align V starting from the left edge, 
    // J starting from the right edge, and align D everywhere.
    //


    /**
     * \brief
     */
    ///@{
    struct CDR3AlignerFunctor_V {
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec, 
                        VDJAlignerParameters params = VDJAlignerParameters()) const
        {
            seq_len_t p_size = pattern.size(), t_size = text.size();
            NoGapAlignment::events_storage_t vec;
            vec.reserve(min(p_size, t_size) + 1);

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                vec.push_back(pattern[i] != text[i]);
            }

            avec->addAlignment(gene, 1, 1, vec);
        }
    };

    struct CDR3AlignerFunctor_D {
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec, 
                        VDJAlignerParameters params = VDJAlignerParameters()) const 
        {
            // seq_len_t match_min_len = params.min_D_len;
            // seq_len_t t_size = text.size(), p_size = pattern.size(), min_size = min(t_size, p_size), min_subsize;
            // seq_len_t p_start, t_start;
            // AlignmentVectorBase::events_storage_t bitvec;
            // bitvec.reserve(t_size + p_size*2);

            // bitvec.push_back(pattern_i[i] == text[text_i + i]);

            // for (seq_len_t text_i = 0; text_i < t_size - match_min_len + 1; ++text_i) {
            //     min_subsize = min((int) p_size, t_size - text_i);
            //     for (seq_len_t i = 0; i < min_subsize; ++i) {
            //         bitvec.push_back(pattern_i[i] == text[text_i + i]);
            //     }
            //     if (open_match && (min_subsize - p_start) >= match_min_len) {
            //         avec->addAlignment(gene, p_start + 1, t_start + 1, min_subsize - p_start);
            //     }
            // }

            // for (seq_len_t pattern_i = 1 /* WTF?! */ ; pattern_i < p_size - match_min_len + 1; ++pattern_i) {
            //     min_subsize = min(p_size - pattern_i, (int) t_size);
            //     for (seq_len_t i = 0; i < min_subsize; ++i) {
            //         if (pattern[pattern_i + i] == text[i]) {
            //             if (!open_match) {
            //                 p_start = pattern_i + i;
            //                 t_start = i;
            //                 open_match = true;
            //             }
            //         } else if (open_match) {
            //             if ((pattern_i + i - p_start) >= match_min_len) {
            //                 avec->addAlignment(gene, p_start + 1, t_start + 1, pattern_i + i - p_start);
            //             }
            //             open_match = false;
            //         }
            //     }
            //     if (open_match && (pattern_i + min_subsize - p_start) >= match_min_len) {
            //         avec->addAlignment(gene, p_start + 1, t_start + 1, pattern_i + min_subsize - p_start);
            //     }
            // }
        }
    };

    struct CDR3AlignerFunctor_J {
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec, 
                        VDJAlignerParameters params = VDJAlignerParameters()) const 
        {
            seq_len_t p_size = pattern.size(), t_size = text.size();
            NoGapAlignment::events_storage_t vec;
            vec.reserve(min(p_size, t_size) + 1);

            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                vec.push_back(pattern[p_size - i - 1] != text[t_size - i - 1]);
            }

            avec->addAlignment(gene, p_size - min(t_size, p_size) + 1, t_size - min(t_size, p_size) + 1, vec);
        }
    };
    ///@}


    /**
     *
     */
    typedef VDJAlignerBase<NoGapAlignmentVector, 
                           CDR3AlignerFunctor_V, 
                           CDR3AlignerFunctor_D, 
                           CDR3AlignerFunctor_J> CDR3NucleotideAligner;


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
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec, 
                        VDJAlignerParameters params = VDJAlignerParameters()) const 
        {
            // avec->addAlignment(gene, );
        }
    };

    struct SWAlignerFunctor_D {
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec, 
                        VDJAlignerParameters params = VDJAlignerParameters()) const 
        {
            // avec->addAlignment(gene, );
        }
    };
    ///@}


    /**
     * \typedef SmithWatermanAligner
     *
     * \brief Smith-Waterman aligner for finding maximal matches with gaps. Don't takes into account errors.
     */
    typedef VDJAlignerBase<GappedAlignmentVector, 
                           SWAlignerFunctor_VJ, 
                           SWAlignerFunctor_D, 
                           SWAlignerFunctor_VJ> SmithWatermanAligner;


    //
    // Smith-Waterman with no gap allowed, but with errors
    //


    /**
     * \brief
     */
    ///@{
    struct SWNGAlignerFunctor_VJ {
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec, 
                        VDJAlignerParameters params = VDJAlignerParameters()) const 
        {
            // avec->addAlignment(gene, );
        }
    };

    struct SWNGAlignerFunctor_D {
        void operator()(seg_index_t gene, 
                        const sequence_t &pattern, 
                        const sequence_t &text, 
                        NoGapAlignmentVector *avec, 
                        VDJAlignerParameters params = VDJAlignerParameters()) const 
        {
            // avec->addAlignment(gene, );
        }
    };
    ///@}


    /**
     * \typedef SmithWatermanNoGapAligner
     *
     * \brief Smith-Waterman aligner without gaps, returns maximal matches with information about mismatch errors.
     */
    typedef VDJAlignerBase<NoGapAlignmentVector,
                           SWNGAlignerFunctor_VJ, 
                           SWNGAlignerFunctor_D, 
                           SWNGAlignerFunctor_VJ> SmithWatermanNoGapAligner;

}

#endif