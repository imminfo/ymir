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


#include <string>
#include <vector>

#include "types.h"
#include "genesegment.h"


using namespace std;


namespace ymir {

    #define DEFAULT_LOCAL_ALIGNMENT_RESERVE 160


    /**
     * \class AbstractVDJAligner
     */
    class AbstractVDJAligner {
    public:

        /**
         * \typedef alignment_score_t
         *
         * \brief Type for scores of gene segments alignments to input sequences.
         */
        typedef int16_t alignment_score_t;


        /**
         * \struct Alignment
         *
         */
        struct Alignment {

            seg_index_t segment;
            seq_len_t start, end;
            alignment_score_t score;
            std::vector<AlignmentEvent> events;


            Alignment()
                    : segment(0), start(0), end(0), score(0)
            { }


            Alignment(seg_index_t segment_, seq_len_t start_, seq_len_t end_, alignment_score_t score_, const std::vector<AlignmentEvent>& events_)
                    : segment(segment_), start(start_), end(end_), score(score_), events(events_)
            { }

        };


        /**
         * \struct GeneSegmentsAlignment
         */
        struct GeneSegmentsAlignment {

            GeneSegmentsAlignment()
                    : _gene(UNDEF_GENE), _n_alignments(0), _alignments(nullptr)
            { }


            GeneSegmentsAlignment(GeneSegments gene, seg_index_t n_alignments, Alignment *alignments)
                    : _gene(gene), _n_alignments(n_alignments), _alignments(alignments)
            { }


            virtual ~GeneSegmentsAlignment() {
//                delete [] _alignments;
//                delete [] _n_alignments;
            }


            GeneSegments gene() const { return _gene; }


            size_t size() const { return _n_alignments; }


            const Alignment& getAlignment(seg_index_t i) const { return _alignments[i]; }

        protected:

            GeneSegments _gene;
            Alignment *_alignments;
            seg_index_t _n_alignments;

        };


        /**
         *
         */
        AbstractVDJAligner(const VDJRecombinationGenes &genes, alignment_score_t threshold) : _threshold(threshold) {
            _genes = new VDJRecombinationGenes(genes);
        }


        virtual ~AbstractVDJAligner() {
            if (_genes) { delete _genes; }
        }


        /**
         * \brief Align the given sequence to the specific gene segment of the specific gene.
         *
         * \param sequence Pattern sequence.
         * \param seg_index Index of the target gene segment.
         */
        ///@{
        virtual Alignment alignVar(const std::string& sequence, seg_index_t seg_index) = 0;

        virtual Alignment alignDiv(const std::string& sequence, seg_index_t seg_index) = 0;

        virtual Alignment alignJoi(const std::string& sequence, seg_index_t seg_index) = 0;
        ///@}


        /**
         * \brief Align the given sequence to all gene segments of the specific gene.
         *
         * \param sequence Pattern sequence.
         */
        ///@{
        virtual void alignVar(const std::string& sequence) {
            std::vector<Alignment> vec;
            vec.reserve(_genes->V().size() / 2 + 2);
            Alignment tmp;
            for (seg_index_t i = 1; i <= _genes->V().max(); ++i) {
//                tmp = this->alignVar(sequence);
//                if (tmp.score >= _threshold) {
//                    vec.push_back(tmp);
//                }
            }

            // ???

        }

        virtual void alignDiv(const std::string &sequence) = 0;

        virtual void alignJoi(const std::string &sequence) = 0;
        ///@}


        /**
         * \brief Access the latest alignment results.
         */
        ///@{
        const GeneSegmentsAlignment& lastVarAlignment() const { return _last_v_alignment; }

        const GeneSegmentsAlignment& lastDivAlignment() const { return _last_d_alignment; }

        const GeneSegmentsAlignment& lastJoiAlignment() const { return _last_j_alignment; }
        ///@}

    protected:

        alignment_score_t _threshold;
        VDJRecombinationGenes *_genes;
        GeneSegmentsAlignment _last_v_alignment, _last_d_alignment, _last_j_alignment;


        AbstractVDJAligner() {
            _genes = nullptr;
        }

    };


//    template <class _Input, class _Output>
    class AbstractAligner {
    public:

        /** \brief Vector of starts and ends of alignment results.
        *
        */
        struct LocalAlignmentIndices {
            // vector of 3-tuples: pattern start, text start, alignment length;
//            _Output *alignment;
            seq_len_t *alignment;
            size_t n;


//            LocalAlignmentIndices(_Output *alignment_, size_t n_) : n(n_ / 4) {
            LocalAlignmentIndices(const vector<seq_len_t> vec) {
                if (vec.size()) {
                    this->n = vec.size() / 3;
                    this->alignment = new seq_len_t[this->n * 3];
                    for (int i = 0; i < this->n * 3; ++i) {
                        this->alignment[i] = vec[i];
                    }
                } else {
                    this->n = 0;
                    this->alignment = nullptr;
                }
            }


            LocalAlignmentIndices(const LocalAlignmentIndices& other) {
                this->n = other.n;
                this->alignment = new seq_len_t[this->n * 3];
                for (int i = 0; i < this->n * 3; ++i) {
                    this->alignment[i] = other.alignment[i];
                }
            }


            virtual ~LocalAlignmentIndices() {
                delete [] this->alignment;
            }


            size_t array_size() const { return this->n * 3; }


            size_t size() const { return this->n; }


            Alignment operator[](size_t index) const {
                return Alignment(this->alignment + 3*index);
            }
        };


//        virtual static const seq_len_t& getEdgeAlignmentMaxIndex(const _Input& pattern, const _Input& text, bool on_reverse_text = false) const =0;


        virtual seq_len_t align5end(const string& pattern, const string& text) const = 0;

        virtual seq_len_t align3end(const string& pattern, const string& text) const = 0;

        virtual LocalAlignmentIndices alignLocal(const string& pattern, const string& text, seq_len_t match_min_len = 3) const = 0;

    };


//    class NaiveNucleotideAligner : public AbstractAligner<string, seq_len_t> {
//    class NaiveNucleotideAligner : public AbstractVDJAligner {
    class NaiveNucleotideAligner : public AbstractAligner {
    public:

        virtual seq_len_t align5end(const string& pattern, const string& text) const {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0;
            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[i] != text[i]) { break; }
                matches += 1;
            }
            return matches;
        }


        virtual seq_len_t align3end(const string& pattern, const string& text) const {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0;
            for (seq_len_t i = 0; i < min(p_size, t_size); ++i) {
                if (pattern[p_size - i - 1] != text[t_size - i - 1]) { break; }
                matches += 1;
            }
            return matches;
        }


        virtual LocalAlignmentIndices alignLocal(const string& pattern, const string& text, seq_len_t match_min_len = 3) const {
            vector<seq_len_t> vec;
            vec.reserve(DEFAULT_LOCAL_ALIGNMENT_RESERVE);

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
//                            vec.push_back(p_start);
//                            vec.push_back(pattern_i + i - 1);
//                            vec.push_back(t_start);
//                            vec.push_back(i - 1);
                            vec.push_back(p_start + 1);
                            vec.push_back(t_start + 1);
                            vec.push_back(pattern_i + i - 1);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (pattern_i + min_subsize - p_start) >= match_min_len) {
//                    vec.push_back(p_start);
//                    vec.push_back(pattern_i + min_subsize - 1);
//                    vec.push_back(t_start);
//                    vec.push_back(min_subsize - 1);
                    vec.push_back(p_start + 1);
                    vec.push_back(t_start + 1);
                    vec.push_back(pattern_i + min_subsize - p_start);
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
//                            vec.push_back(p_start);
//                            vec.push_back(i - 1);
//                            vec.push_back(t_start);
//                            vec.push_back(text_i + i - 1);
                            vec.push_back(p_start + 1);
                            vec.push_back(t_start + 1);
                            vec.push_back(i - 1);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (min_subsize - p_start) >= match_min_len) {
//                    vec.push_back(p_start);
//                    vec.push_back(min_subsize - 1);
//                    vec.push_back(t_start);
//                    vec.push_back(text_i + min_subsize - 1);
                    vec.push_back(p_start + 1);
                    vec.push_back(t_start + 1);
                    vec.push_back(min_subsize - p_start);
                }
            }

//            for (size_t i = 0; i < vec.size(); ++i) { ++vec[i]; }
            return LocalAlignmentIndices(vec);
        }

    };


//    class NaiveAminoAcidAligner : public AbstractAligner<string, seq_len_t> {
    class NaiveAminoAcidAligner : public AbstractAligner {

    public:

        NaiveAminoAcidAligner() { }


        virtual seq_len_t align5end(const string& pattern, const string& text) const {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0, max_matches = 0, all_matches = 0;
            string codon_s = "";
            for (seq_len_t i = 0; i < std::min((seq_len_t) (1 + p_size / 3), t_size); ++i) {
                max_matches = 0;
                // go through all codons and find the maximal match
                CodonTable::Codons codon = _codons.codons(text[i]);
                while(!codon.end()) {
                    matches = 0;
                    codon_s = codon.next();
                    if (pattern[i*3] == codon_s[0] && i*3 < p_size) {
                        ++matches;
                        if (pattern[i*3 + 1] == codon_s[1] && (i*3 + 1) < p_size) {
                            ++matches;
                            if (pattern[i*3 + 2] == codon_s[2] && (i*3 + 2) < p_size) {
                                ++matches;
                            }
                        }
                    }

                    if (matches > max_matches) {
                        max_matches = matches;
                        if (max_matches == 3) { break; }
                    }
                }

                // if match == 3 then go to the next amino acid
                all_matches += max_matches;
                if (max_matches != 3) { break; }
            }

            return all_matches;
        }


        virtual seq_len_t align3end(const string& pattern, const string& text) const {
            seq_len_t p_size = pattern.size(), t_size = text.size(), matches = 0, max_matches = 0, all_matches = 0;
            string codon_s = "";
            for (seq_len_t i = 0; i < std::min((seq_len_t) (1 + p_size / 3), t_size); ++i) {
                max_matches = 0;
                // go through all codons and find the maximal match
                CodonTable::Codons codon = _codons.codons(text[t_size - i - 1]);
                while(!codon.end()) {
                    matches = 0;
                    codon_s = codon.next();
                    if (pattern[p_size - 1 - i*3] == codon_s[2] && i*3 < p_size) {
                        ++matches;
                        if (pattern[p_size - 1 - i*3 - 1] == codon_s[1] && (i*3 + 1) < p_size) {
                            ++matches;
                            if (pattern[p_size - 1 - i*3 - 2] == codon_s[0] && (i*3 + 2) < p_size) {
                                ++matches;
                            }
                        }
                    }

                    if (matches > max_matches) {
                        max_matches = matches;
                        if (max_matches == 3) { break; }
                    }
                }

                // if match == 3 then go to the next amino acid
                all_matches += max_matches;
                if (max_matches != 3) { break; }
            }

            return all_matches;
        }


        virtual LocalAlignmentIndices alignLocal(const string& pattern, const string& text, seq_len_t match_min_len = 3) const {
            return LocalAlignmentIndices(vector<seq_len_t>(1));
        }

    protected:

        const CodonTable _codons;

    };


//    class NoGapNucleotideAligner : public AbstractVDJAligner {
//
//    };
}

#endif