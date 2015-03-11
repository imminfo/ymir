/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vadim dot nazarov at mailbox dot com>
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

#ifndef _MAAG_H
#define _MAAG_H


#include "multimatrixchain.h"


namespace ymir {


    /**
    * \class MAAG
    *
    * \brief Multi-Alignment Assembly Graph - basic class for representing
    * all possible generation scenario of a nucleotide or amino acid sequence of an immune receptor.
    */
    class MAAG : protected ProbMMC {

    public:

        MAAG() {
            _events = nullptr;
        }


        MAAG(bool fill_event_indices) : MAAG() {
            if (fill_event_indices) {
                _events = new EventIndMMC();
            }
        }


        virtual ~MAAG() {
            if (_events) { delete _events; }
        }


        /**
        * \brief Compute and return full assembling probability of this sequence with specific gene segments alignments.
        *
        * \return Slice assembling probability.
        */
        numeric sliceProbability() const {
//            matrix_t res(this->_chain[0]);
//            for (matrix_ind i = 1; i < this->_chain.size(); i++) {
//                res = res * this->_chain[i];
//            }
//            return res(0,0);
        }


        /**
        * \brief Compute and return full assembling probability of this sequence (i.e., with all gene segments alignments).
        *
        * \return Full assembling probability.
        */
        numeric fullProbability() const {
//            matrix_t res(this->_chain[0]);
//            for (matrix_ind i = 1; i < this->_chain.size(); i++) {
//                res = res * this->_chain[i];
//            }
//            return res(0,0);
        }


    protected:

        EventIndMMC *_events;  /** Matrix of indices of events for each edge. */

        // gene segment metadata structure
        eventind_t *_vind, *_jind, *_dind;
        prob_t *_vprobs, *_jprobs, *_dprobs;  // ??? DJ joint probability;
        eventind_t _n_v, _n_j, _n_d;

        seq_len_t *_seq_poses;  /** Vector of the initial clonal sequence's positions for each vertex. */
        seq_len_t _seq_len;  /** Size of the initial clonal sequence. */

    };
}

#endif