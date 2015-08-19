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
#include "math.h"
#include <iomanip>

#include "types.h"
//#include "repertoire.h"


namespace ymir {


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


    std::string translate(const std::string &nuc_seq) {
        return "nothing";
    }


    prob_t loglikelihood(const std::vector<prob_t> &vec, prob_t laplace = 1e-80) {
        prob_t res = 0;
        for (size_t i = 0; i < vec.size(); ++i) { res += log10((vec[i] + laplace)); } // what to do in case of high precision numbers?
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
            if (loglikelihood(prob_vec) > prev_ll) { std::cout << " (grows)"; }
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
        std::string tmp;
        for (size_t i = 0; i < sequence.size(); i += 3) {
            tmp = sequence.substr(i, 3);
            if (tmp == "TAG" || tmp == "TAA" || tmp == "TGA") {
                return true;
            }
        }
        return false;
    }

}

#endif