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
    public:

        /**
         *
         */
        MAAG(const MAAG &other) {
            if (other._events) {
                _events = new EventIndMMC(*other._events);
            }
            if (other._sequence) { _sequence = new string(*other._sequence); }
        }


        /**
         * \brief Special swap constructor for MAAGs that will be used only for computation of
         * the full probability.
         */
        MAAG(ProbMMC &prob_mcc) {
            this->swap(prob_mcc);
            _events = nullptr;
            _sequence = nullptr;
            _seq_poses = nullptr;
            _n_poses = 0;
        }


        /**
         * \brief Special swap constructor for MAAGs that will be used for statistical inference.
         */
        MAAG(ProbMMC &prob_mcc, EventIndMMC &eventind_mcc, string sequence, seq_len_t *seq_poses, seq_len_t n_poses) {
            this->swap(prob_mcc);
            _events = new EventIndMMC();
            _events->swap(eventind_mcc);
            _sequence = new string(sequence);
            _seq_poses = seq_poses;
            _n_poses = n_poses;
        }


        virtual ~MAAG() {
            if (_events) { delete _events; }
            if (_seq_poses) { delete [] _seq_poses; }
            if (_sequence) { delete _sequence; }
        }


        /**
        * \brief Compute and return full assembling probability of this sequence (i.e., with all gene segments alignments).
        *
        * \return Full assembling probability.
        */
        ///@{
        numeric fullProbability(eventind_t v_index, eventind_t j_index) const {
            // P(Vi) * P(#dels | Vi) * P(V-J insertion seq) * P(#dels | Ji) * P(Ji)
            return (_chain[0][v_index] *        // P(Vi)
                    _chain[1][v_index] *        // P(#dels | Vi)
                    _chain[2][0] *              // P(V-J insertion seq)
                    _chain[3][j_index] *        // P(#dels | Ji)
                    _chain[4][j_index])(0, 0);  // P(Ji)
        }
        numeric fullProbability(eventind_t v_index, eventind_t d_index, eventind_t j_index) const {
            // P(Vi) * P(#dels | Vi) * P(V-D3' insertion seq) * P(D5'-D3' deletions | Di) * P(D5'-J insertion seq) * P(#dels | Ji) * P(Ji & Di)
            return (_chain[0][v_index] *      // P(Vi)
                    _chain[1][v_index] *      // P(#dels | Vi)
                    _chain[2][0] *            // P(V-D3' insertion seq)
//                    _chain[2][d_index] *            // P(V-D3' insertion seq)
                    _chain[3][d_index] *      // P(D5'-D3' deletions | Di)
                    _chain[4][0] *            // P(D5'-J insertion seq)
//                    _chain[4][d_index] *            // P(D5'-J insertion seq)
                    _chain[5][j_index] *      // P(#dels | Ji)
                    _chain[6][0](j_index, d_index))(0, 0);  // P(Ji & Di)
//                    _jdata->prob(j_index * _ddata->size() + d_index))(0, 0);  // P(Ji & Di)
        }
        ///@}


        ///@{
        eventind_t nVar() { return _chain[0].size(); }
        eventind_t nJoi() { return _chain[_chain.size() - 2].size(); }
        eventind_t nDiv() { return (_chain.size() == 5) ? 0 : _chain[3].size(); }
        ///@}


        ///@{
//        segindex_t getVar(uint8_t i) { return _events ? (*_events)[0][i] : 0; }
//        segindex_t getJoi(uint8_t i) { return _events ? (*_events)[_chain.size() - 2][i] : 0; }
//        segindex_t getDiv(uint8_t i) { return 0; } // ???
        ///@}


        /**
         * \brief Access to event indices and event probabilities in the underlying matrices of the MAAG.
         *
         * \param node_i Index of the node.
         * \param mat_i Index of the matrix in the node.
         * \param row Which row to choose.
         * \param col Which column to choose.
         *
         * \return Event probability or event index. In the second case zero will be returned if no event chain matrix
         * is stored in this MAAG.
         */
        ///@{
        prob_t event_probability(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
            return (*this)(node_i, mat_i, row, col);
        }
        eventind_t event_index(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
            return _events ? (*_events)(node_i, mat_i, row, col) : 0;
        };
        ///@}


        seq_len_t position(seq_len_t i) const { return _seq_poses[i]; }


        void serialize(const string& filepath) const {
            
        }


    protected:

        EventIndMMC *_events;  /** Matrix of indices of events for each edge. */

        seq_len_t *_seq_poses;  /** Vector of the initial clonotype sequence's positions for each vertex. */
        seq_len_t _n_poses;
        string *_sequence;  /** Nucleotide or amino acid CDR3 sequence. */


        MAAG() {
            _events = nullptr;
            _sequence = nullptr;
            _seq_poses = nullptr;
        }

    };
}

#endif