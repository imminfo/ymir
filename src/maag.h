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

#ifndef _MAAG_H_
#define _MAAG_H_


#include "multimatrixchain.h"


namespace ymir {

    class MAAG;


    /**
    * \class MAAG
    *
    * \brief Multi-Alignment Assembly Graph - basic class for representing
    * all possible generation scenarios of a nucleotide or amino acid sequence of an immune receptor.
    */
    class MAAG : protected ProbMMC {

    protected:

        /**
        * \struct GeneSegment
        */
        struct GeneSegment {

        public:

            GeneSegment(eventind_t n) {
                _n = n;
                _indices = new eventind_t[_n];
            }


            virtual ~GeneSegment() {
                delete [] _indices;
            }


            eventind_t event_index(eventind_t i) const { return _indices[i]; }
            eventind_t& event_index(eventind_t i) { return _indices[i]; }


            eventind_t size() const { return _n; }

        private:

            eventind_t _n;
            eventind_t *_indices;

        };

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
            // delete [] _seq_poses;
        }


        /**
        * \brief Compute and return full assembling probability of this sequence (i.e., with all gene segments alignments).
        *
        * \return Full assembling probability.
        */
        numeric fullProbability(eventind_t v_index, eventind_t j_index = 0) const {
            // P(Vi) * P(#dels | Vi) * P(V-J insertion seq) * P(#dels | Ji) * P(Ji)
            return (_chain[0][v_index] *        // P(Vi)
                    _chain[1][v_index] *        // P(#dels | Vi)
                    _chain[2][0] *              // P(V-J insertion seq)
                    _chain[3][j_index] *        // P(#dels | Ji)
                    _chain[4][j_index])(0, 0);  // P(Ji)
        }
        numeric fullProbability(eventind_t v_index, eventind_t d_index, eventind_t j_index = 0) const {
            // P(Vi) * P(#dels | Vi) * P(V-D3' insertion seq) * P(D5'-D3' deletions | Di) * P(D5'-J insertion seq) * P(#dels | Ji) * P(Ji & Di)
            return (_chain[0][v_index] *      // P(Vi)
                    _chain[1][v_index] *      // P(#dels | Vi)
                    _chain[2][0] *            // P(V-D3' insertion seq)
                    _chain[3][d_index] *      // P(D5'-D3' deletions | Di)
                    _chain[4][0] *            // P(D5'-J insertion seq)
                    _chain[5][j_index] *      // P(#dels | Ji)
                    _chain[6][0](j_index, d_index))(0, 0);  // P(Ji & Di)
//                    _jdata->prob(j_index * _ddata->size() + d_index))(0, 0);  // P(Ji & Di)
        }


        ///@{
        eventind_t nV() { return _vdata->size(); }
        eventind_t nJ() { return _jdata->size(); }
        eventind_t nD() { return _ddata->size(); }
        ///@}


        ///@{
        eventind_t vgene(uint8_t i) { return _vdata->event_index(i); }
        eventind_t jgene(uint8_t i) { return _jdata->event_index(i); }
        eventind_t dgene(uint8_t i) { return _ddata->event_index(i); }
        ///@}



    protected:

        EventIndMMC *_events;  /** Matrix of indices of events for each edge. */

        GeneSegment *_vdata, *_jdata, *_ddata;

        seq_len_t *_seq_poses;  /** Vector of the initial clonal sequence's positions for each vertex. */
        seq_len_t _seq_len;  /** Size of the initial clonal sequence. */

    };
}

#endif