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


using namespace std;


namespace ymir {

    #define DEFAULT_LOCAL_ALIGNMENT_RESERVE 160


//    template <class _Input, class _Output>
    class AbstractAligner {
    public:

        /** \brief Vector of starts and ends of alignment results.
        *
        */
        struct LocalAlignmentIndices {
            // vector of 4-tuples: pattern start, pattern end, text start, text end;
//            _Output *alignment;
            seq_len_t *alignment;
            size_t n;


//            LocalAlignmentIndices(_Output *alignment_, size_t n_) : n(n_ / 4) {
            LocalAlignmentIndices(const vector<seq_len_t> vec) {
                if (vec.size()) {
                    this->n = vec.size() / 4;
                    this->alignment = new seq_len_t[this->n * 4];
                    for (int i = 0; i < this->n * 4; ++i) {
                        this->alignment[i] = vec[i];
                    }
                } else {
                    this->n = 0;
                    this->alignment = nullptr;
                }
            }


            LocalAlignmentIndices(const LocalAlignmentIndices& other) {
                this->n = other.n;
                this->alignment = new seq_len_t[this->n * 4];
                for (int i = 0; i < this->n * 4; ++i) {
                    this->alignment[i] = other.alignment[i];
                }
            }


            virtual ~LocalAlignmentIndices() {
                delete [] this->alignment;
            }


            size_t array_size() const { return this->n * 4; }


            size_t size() const { return this->n; }


            d_alignment_t operator[](size_t index) const {
                return d_alignment_t(this->alignment + 4*index);
            }
        };


//        virtual static const seq_len_t& getEdgeAlignmentMaxIndex(const _Input& pattern, const _Input& text, bool on_reverse_text = false) const =0;


        virtual seq_len_t align5end(const string& pattern, const string& text) const = 0;

        virtual seq_len_t align3end(const string& pattern, const string& text) const = 0;

        virtual LocalAlignmentIndices alignLocal(const string& pattern, const string& text, seq_len_t match_min_len = 3) const = 0;

    };


//    class NaiveNucleotideAligner : public AbstractAligner<string, seq_len_t> {
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
                            vec.push_back(p_start);
                            vec.push_back(pattern_i + i - 1);
                            vec.push_back(t_start);
                            vec.push_back(i - 1);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (pattern_i + min_subsize - p_start) >= match_min_len) {
                    vec.push_back(p_start);
                    vec.push_back(pattern_i + min_subsize - 1);
                    vec.push_back(t_start);
                    vec.push_back(min_subsize - 1);
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
                            vec.push_back(p_start);
                            vec.push_back(i - 1);
                            vec.push_back(t_start);
                            vec.push_back(text_i + i - 1);
                        }
                        open_match = false;
                    }
                }
                if (open_match && (min_subsize - p_start) >= match_min_len) {
                    vec.push_back(p_start);
                    vec.push_back(min_subsize - 1);
                    vec.push_back(t_start);
                    vec.push_back(text_i + min_subsize - 1);
                }
            }

            for (size_t i = 0; i < vec.size(); ++i) { ++vec[i]; }
            return LocalAlignmentIndices(vec);
        }

    };


//    class NaiveAminoAcidAligner : public AbstractAligner<string, seq_len_t> {
    class NaiveAminoAcidAligner : public AbstractAligner {
    public:

        virtual seq_len_t align5end(const string& pattern, const string& text) const {
            return -1;
        }


        virtual seq_len_t align3end(const string& pattern, const string& text) const {
            return -1;
        }


        virtual LocalAlignmentIndices alignLocal(const string& pattern, const string& text, seq_len_t match_min_len = 3) const {
            return LocalAlignmentIndices(vector<seq_len_t>(1));
        }

    };
}

#endif