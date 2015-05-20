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

#define VJ_CHAIN_SIZE 4
#define VDJ_CHAIN_SIZE 7


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

        friend class MAAGForwardBackwardAlgorithm;  // ):

    public:

        /**
         *
         */
        MAAG(const MAAG &other) : ProbMMC(other) {
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


        /**
         * 
         */
        virtual ~MAAG() {
            if (_events) { delete _events; }
            if (_seq_poses) { delete [] _seq_poses; }
            if (_sequence) { delete _sequence; }
        }


        /**
         * \brief Compute and return full assembling probability of this sequence (i.e., with all gene segments alignments).
         *
         * \param v_index Which aligned Variable gene segment to choose.
         * \param d_index Which aligned Diversity gene segment to choose.
         * \param j_index Which aligned Joining gene segment to choose.
         * \param action
         *
         * \return Full assembling probability. Returns 0 if wrong recombination function is used.
         */
        ///@{
        prob_t fullProbability(eventind_t v_index, eventind_t j_index) const {
            if (this->is_vj()) {
                // P(Vi) * P(#dels | Vi) * P(V-J insertion seq) * P(#dels | Ji) * P(Ji)
                return (_chain[0][0](v_index, j_index) *  // P(Vi & Ji)
                        _chain[1][v_index] *        // P(#dels | Vi)
                        _chain[2][0] *              // P(V-J insertion seq)
                        _chain[3][j_index])         // P(#dels | Ji)
                        .sum();
            } else {
                return 0;
            }
        }

        prob_t fullProbability(eventind_t v_index, eventind_t d_index, eventind_t j_index) const {
            if (this->is_vdj()) {
                // P(Vi) * P(#dels | Vi) * P(V-D3' insertion seq) * P(D5'-D3' deletions | Di) * P(D5'-J insertion seq) * P(#dels | Ji) * P(Ji & Di)
                return (_chain[0][v_index] *      // P(Vi)
                        _chain[1][v_index] *      // P(#dels | Vi)
                        _chain[2][0] *            // P(V-D3' insertion seq)
                        _chain[3][d_index] *      // P(D5'-D3' deletions | Di)
                        _chain[4][0] *            // P(D5'-J insertion seq)
                        _chain[5][j_index] *      // P(#dels | Ji)
                        _chain[6][0](j_index, d_index))(0, 0);  // P(Ji & Di)
            } else {
                return 0;
            }
        }

        prob_t fullProbability(MAAG_COMPUTE_PROB_ACTION action = SUM_PROBABILITY) const {
            // choose the max full probability from all possible recombinations of V(D)J gene segment indices
            if (action == MAX_PROBABILITY) {
                prob_t max_prob = 0, cur_prob = 0;
                if (this->is_vj()) {
                    for (eventind_t v_index = 0; v_index < this->nVar(); ++v_index) {
                        for (eventind_t j_index = 0; j_index < this->nJoi(); ++j_index) {
                            cur_prob = this->fullProbability(v_index, j_index);
                            if (cur_prob > max_prob) { max_prob = cur_prob; }
                        }
                    }
                } else {
                    for (eventind_t v_index = 0; v_index < this->nVar(); ++v_index) {
                        for (eventind_t d_index = 0; d_index < this->nDiv(); ++d_index) {
                            for (eventind_t j_index = 0; j_index < this->nJoi(); ++j_index) {
                                cur_prob = this->fullProbability(v_index, d_index, j_index);
                                if (cur_prob > max_prob) { max_prob = cur_prob; }
                            }
                        }
                    }
                }
                return max_prob;
            }
            // compute the sum of full probabilities of all possible recombinations of V(D)J gene segment indices
            else {
                prob_t sum_prob = 0;
                if (this->is_vj()) {
                    for (eventind_t v_index = 0; v_index < this->nVar(); ++v_index) {
                        for (eventind_t j_index = 0; j_index < this->nJoi(); ++j_index) {
                            sum_prob += this->fullProbability(v_index, j_index);
                        }
                    }
                } else {
                    for (eventind_t v_index = 0; v_index < this->nVar(); ++v_index) {
                        for (eventind_t d_index = 0; d_index < this->nDiv(); ++d_index) {
                            for (eventind_t j_index = 0; j_index < this->nJoi(); ++j_index) {
                                sum_prob += this->fullProbability(v_index, d_index, j_index);
                            }
                        }
                    }
                }
                return sum_prob;
            }
        }
        ///@}


        /**
         * \brief Get the number of aligned gene segments.
         *
         * \return Number of the aligned specific gene segments.
         */
        ///@{
        eventind_t nVar() const { return _chain[0].size(); }
        eventind_t nJoi() const { return _chain[_chain.size() - 2].size(); }
        eventind_t nDiv() const { return (_chain.size() == VJ_CHAIN_SIZE) ? 0 : _chain[3].size(); }
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


        /**
         * \brief Save the graph to the harddrive and load a previously saved graph.
         *
         * \param filepath Path to the file with a saved graph for loading and new file for saving.
         * \param stream An output / input stream, from that read the graph or save the graph.
         */
        ///@{
        bool save(const string &filepath) const {
            ofstream ofs(filepath);
            return this->save(ofs);
        }

        bool save(ostream &stream) const {
            return false;
        }

        bool load(const string &filepath) {
            ifstream ifs(filepath);
            return this->load(ifs);
        }

        bool load(istream &stream) {
            return false;
        }
        ///@}


        /**
         * \brief Get is this MAAG was built from clonotype either with VJ or with VDJ recombination.
         */
        ///@{
        bool is_vj() const { return _chain.size() == VJ_CHAIN_SIZE; }
        bool is_vdj() const { return _chain.size() == VDJ_CHAIN_SIZE; }
        ///@}

        RECOMBINATION recombination() const { return _chain.size() == VJ_CHAIN_SIZE ? VJ_RECOMB : VDJ_RECOMB ; }

    protected:

        EventIndMMC *_events;  /** Matrix of indices of events for each edge. */

        seq_len_t *_seq_poses;  /** Vector of the initial clonotype sequence's positions for each vertex. */
        seq_len_t _n_poses;
        string *_sequence;  /** Nucleotide or amino acid CDR3 sequence. */



        /**
         * \brief Default constructor.
         */
        MAAG() {
            _events = nullptr;
            _sequence = nullptr;
            _seq_poses = nullptr;
        }

    };
}

#endif