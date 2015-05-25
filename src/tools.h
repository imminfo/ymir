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


    prob_t loglikelihood(const std::vector<prob_t> &vec) {
        prob_t res = 0, vec_sum = 0;
        for (size_t i = 0; i < vec.size(); ++i) { vec_sum += vec[i]; }
        for (size_t i = 0; i < vec.size(); ++i) { res += std::log10(vec[i] / vec_sum); } // what to do in case of high precision numbers?
        return res;
    }

}

#endif