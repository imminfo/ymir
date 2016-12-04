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

#ifndef _TOOLS_H_
#define _TOOLS_H_


#include <string>
#include <vector>
#include <chrono>
#include <ctime>
#include "math.h"

#include "types.h"


namespace ymir {


    inline void check_and_throw(bool condition, const std::string &error_msg) {
#ifndef DNDEBUG
        if (condition) {
             throw(std::runtime_error(error_msg));
        }
#endif
    }


    inline uint8_t nuc_hash(char nuc) {
        switch (nuc) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default : return 4;
        }
    }


    inline char inv_nuc_hash(uint8_t hash) {
        switch (hash) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default : return 4;
        }
    }


    inline char complement(char nuc) {
        switch (nuc) {
            case 'A': return 'T';
            case 'C': return 'G';
            case 'G': return 'C';
            case 'T': return 'A';
            default : return NULL_CHAR;
        }
    }


    prob_t loglikelihood(const std::vector<prob_t> &vec, prob_t laplace = 1e-80) {
        prob_t res = 0;
        for (size_t i = 0; i < vec.size(); ++i) {
//            res += log10((vec[i] + laplace));

            if (vec[i] && !std::isnan(vec[i])) {
                res += log10(vec[i]);
            }

//            if ((vec[i] > 1e-60) && !std::isnan(vec[i])) {
//                res += log10(vec[i]);
//            }
        } // what to do in case of high precision numbers?
        // w/ laplace and w/o zero check -77603.267; -100500.399
        // w/o laplace and w/ zero check -23443.267; -46340.3987
        // w/o laplace and w/ eps check  -23443.267; -46340.3987
        return res;
    }


    void prob_summary(const std::vector<prob_t> &prob_vec, prob_t prev_ll = 0) {
        size_t zeros = 0, negative = 0, bignums = 0, nans = 0;

        for (size_t i = 0; i < prob_vec.size(); ++i) {
            zeros += prob_vec[i] ? 0 : 1;
            negative += prob_vec[i] < 0 ? 1 : 0;
            bignums += prob_vec[i] > 1 ? 1 : 0;
            nans += std::isnan(prob_vec[i]) ? 1 : 0;
        }

        std::cout << "Loglikelihood:\t" << std::setprecision(9) << loglikelihood(prob_vec);
        if (prev_ll) {
            if (std::abs(loglikelihood(prob_vec) - prev_ll) < 1e-5) { std::cout << " (exact)"; }
            else if (loglikelihood(prob_vec) > prev_ll) { std::cout << " (grows)"; }
            else { std::cout << " (drops)"; }
        }
        std::cout << std::endl;
        std::cout << "Error probabilities:\t" << (size_t) (zeros + negative + bignums + nans) << std::endl;
        if (zeros) std::cout << "  Zeros:            \t" << (size_t) zeros << std::endl;
        if (nans) std::cout << "  NaNs:             \t" << (size_t) nans << std::endl;
        if (negative) std::cout << "  Negatives:        \t" << (size_t) negative << std::endl;
        if (bignums) std::cout << "  Bigger than one:  \t" << (size_t) bignums << std::endl;
    }


    inline bool is_out_of_frame(const std::string &sequence) {
        return sequence.size() % 3 != 0;
    }


    inline bool has_end_codon(const std::string &sequence) {
        for (size_t i = 0; i < sequence.size(); i += 3) {
            if (sequence[i] == 'T'
                && ((sequence[i+1] == 'A' && sequence[i+2] == 'G')
                    || (sequence[i+1] == 'A' && sequence[i+2] == 'A')
                    || (sequence[i+1] == 'G' && sequence[i+2] == 'A')))
            {
                return true;
            }
        }
        return false;
    }

    inline bool has_bad_aa_codons(const sequence_t &sequence) {
        for (size_t i = 0; i < sequence.size(); ++i) {
            if (sequence[i] == '~' || sequence[i] == '*') {
                return true;
            }
        }
        return false;
    }


    inline bool has_oof_aa_codon(const sequence_t &sequence) {
        for (size_t i = 0; i < sequence.size(); ++i) {
            if (sequence[i] == '~') {
                return true;
            }
        }
        return false;
    }


    inline std::string time_format(std::time_t some_time) {
        struct std::tm *ptm = std::gmtime(&some_time);
        char buffer[30];
        std::strftime(buffer, 30, "%Hh%Mm%Ss", ptm);
//        delete ptm;
        return std::string(buffer);
    }


    inline std::string time_diff(std::chrono::system_clock::time_point tp1, std::chrono::system_clock::time_point tp2) {
        std::time_t diff = std::chrono::system_clock::to_time_t(tp2) - std::chrono::system_clock::to_time_t(tp1);
        return time_format(diff);
    }


    inline std::string time_diff_now(std::chrono::system_clock::time_point tp1) {
        return time_diff(tp1, std::chrono::system_clock::now());
    }


    inline std::string est_time(std::chrono::system_clock::time_point tp1, size_t n_objects, size_t cur_object) {
        return time_format(static_cast<std::time_t>(((n_objects - cur_object) / (float) cur_object
                                                     * (std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) - std::chrono::system_clock::to_time_t(tp1)))));
    }


    inline std::string tot_time(std::chrono::system_clock::time_point tp1, size_t n_objects, size_t cur_object) {
        return time_format(static_cast<std::time_t>(((std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) - std::chrono::system_clock::to_time_t(tp1)) / ((float) cur_object / n_objects))));
    }


    inline std::string print_time(std::chrono::system_clock::time_point tp1, size_t n_objects, size_t cur_object) {
        return std::string("(pass.: ") + time_diff_now(tp1) + std::string("; est.: ") + est_time(tp1, n_objects, cur_object) + std::string("; tot.: ") + tot_time(tp1, n_objects, cur_object) + std::string(")");
    }


    inline seq_len_t pos_codon0(seq_len_t nuc_pos) {
        return (nuc_pos - 1) / 3;
    }

    inline seq_len_t in_codon0(seq_len_t nuc_pos) {
        return (nuc_pos - 1) % 3;
    }

}

#endif