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

#include "types.h"


namespace ymir {


    void write_matrix(const std::string &filepath) {

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


    std::string translate(const std::string &nuc_seq) {
        return "nothing";
    }


    uint sample(uint vec_size, prob_t *probs) {

    }


    prob_t loglikelihood(const std::vector<prob_t> &vec, prob_t laplace = 1e-50) {
        prob_t res = 0, vec_sum = laplace * vec.size();
        for (size_t i = 0; i < vec.size(); ++i) { vec_sum += vec[i]; }
        for (size_t i = 0; i < vec.size(); ++i) { res += std::log10((vec[i] + laplace) / vec_sum); } // what to do in case of high precision numbers?
        return res;
    }


    void prob_summary(const std::vector<prob_t> &prob_vec) {
        size_t zeros = 0, negative = 0, bignums = 0, nans = 0;

        for (size_t i = 0; i < prob_vec.size(); ++i) {
            zeros += prob_vec[i] ? 0 : 1;
            negative += prob_vec[i] < 0 ? 1 : 0;
            bignums += prob_vec[i] > 1 ? 1 : 0;
            nans += isnan(prob_vec[i]) ? 1 : 0;
        }

        std::cout << "Loglikelihood:\t" << loglikelihood(prob_vec) << std::endl;
        std::cout << "Error probabilities:\t" << (size_t) (zeros + negative + bignums + nans) << std::endl;
        std::cout << "  Zeros:\t" << (size_t) zeros << std::endl;
        std::cout << "  NaNs:\t" << (size_t) nans << std::endl;
        std::cout << "  Negatives:\t" << (size_t) negative << std::endl;
        std::cout << "  Bigger than one:\t" << (size_t) bignums << std::endl;
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