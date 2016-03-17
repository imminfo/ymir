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
//    class MAAG : public ProbMMC {
    class MAAG : protected ProbMMC {

        friend class MAAGForwardBackwardAlgorithm;  // ):
        friend class MAAGBuilder;  // )):

    public:

        /**
         * \brief Default constructor.
         */
        MAAG()
            : _recomb(UNDEF_RECOMB),
              _events(nullptr),
              _errors(nullptr),
              _sequence(nullptr),
              _seq_poses(nullptr),
              _n_poses(0),
              _seq_type(NUCLEOTIDE),
              _err_prob(0)
        {
            _values.reserve(1);
        }

        /**
         *
         */
        MAAG(const MAAG &other)
            : ProbMMC(other),
              _recomb(other._recomb),
              _n_poses(other._n_poses),
              _seq_type(other._seq_type),
              _events(other._events ? new EventIndMMC(*other._events) : nullptr),
              _errors(other._errors ? new ErrMMC(*other._errors) : nullptr),
              _err_prob(other._err_prob)
        {
            if (other._seq_poses) {
                _seq_poses.reset(new seq_len_t[_n_poses]);
                std::copy(other._seq_poses.get(), other._seq_poses.get() + _n_poses, _seq_poses.get());
            } else {
                _seq_poses = nullptr;
            }

            if (other._sequence) { _sequence.reset(new std::string(*other._sequence)); }
            else { _sequence = nullptr; }
        }


        MAAG(MAAG &&other) {
            std::swap(_recomb, other._recomb);

            _chain.swap(other._chain);
            _values.swap(other._values);

            _events.swap(other._events);

            _errors.swap(other._errors);
            std::swap(_err_prob, other._err_prob);

            std::swap(_n_poses, other._n_poses);
            _seq_poses.swap(other._seq_poses);

            _sequence.swap(other._sequence);
            std::swap(_seq_type, other._seq_type);
        }


        /**
         * \brief Special swap constructor for MAAGs that will be used only for computation of
         * the full probability.
         */
//        MAAG(ProbMMC &prob_mmc)
//            : _recomb(prob_mmc.chainSize() == VJ_CHAIN_SIZE ? VJ_RECOMB : VDJ_RECOMB),  // TODO: fix this
//              _events(nullptr),
//              _sequence(nullptr),
//              _seq_poses(nullptr),
//              _n_poses(0),
//              _seq_type(NUCLEOTIDE)
//        {
//            this->swap(prob_mmc);
//        }
//
//
//        MAAG(ProbMMC &prob_mmc, ErrMMC &err_mmc)
//            : MAAG(prob_mmc)
//        {
//            _errors.reset(new ErrMMC());
//            _errors.get()->swap(err_mmc);
//        }


        /**
         * \brief Special swap constructor for MAAGs that will be used for statistical inference.
         */
//        MAAG(ProbMMC &prob_mmc,
//             EventIndMMC &eventind_mcc,
//             const std::string &sequence,
//             unique_ptr<seq_len_t[]> &seq_poses,
//             seq_len_t n_poses,
//             SequenceType seq_type)
//                : _recomb(prob_mmc.chainSize() == VJ_CHAIN_SIZE ? VJ_RECOMB : VDJ_RECOMB),
//                  _n_poses(n_poses),
//                  _seq_type(seq_type)
//        {
//            this->swap(prob_mmc);
//
//            _events.reset(new EventIndMMC());
//            _events.get()->swap(eventind_mcc);
//
//            _sequence.reset(new std::string(sequence));
//
//            _seq_poses.swap(seq_poses);
//        }
//
//
//        MAAG(ProbMMC &prob_mmc,
//             EventIndMMC &eventind_mcc,
//             ErrMMC &err_mmc,
//             const std::string &sequence,
//             unique_ptr<seq_len_t[]> &seq_poses,
//             seq_len_t n_poses,
//             SequenceType seq_type)
//            : MAAG(prob_mmc, eventind_mcc, sequence, seq_poses, n_poses, seq_type)
//        {
//            _errors.reset(new ErrMMC());
//            _errors.get()->swap(err_mmc);
//        }


        /**
         *
         */
        virtual ~MAAG() { }


        MAAG& operator= (const MAAG &other) {
            _recomb = other._recomb;
            _chain = other._chain;
            _values = other._values;

            if (other._errors) {
                _errors.reset(new ErrMMC(*other._errors));
            } else {
                _errors.reset();
            }
            _err_prob = other._err_prob;

            if (other._events) {
                _events.reset(new EventIndMMC(*other._events));
            }
            else { _events.reset(); }

            _n_poses = other._n_poses;

            if (other._seq_poses) {
                _seq_poses.reset(new seq_len_t[_n_poses]);
                std::copy(other._seq_poses.get(), other._seq_poses.get() + _n_poses, _seq_poses.get());
            } else {
                _seq_poses = nullptr;
            }

            if (other._sequence) {
                _sequence.reset(new std::string(*other._sequence));
            }
            else { _sequence.reset(); }

            _seq_type = other._seq_type;

            return *this;
        }


        MAAG& operator=(MAAG &&other) {
            std::swap(_recomb, other._recomb);

            _chain.swap(other._chain);
            _values.swap(other._values);

            _events.swap(other._events);

            _errors.swap(other._errors);
            std::swap(_err_prob, other._err_prob);

            std::swap(_n_poses, other._n_poses);
            _seq_poses.swap(other._seq_poses);

            _sequence.swap(other._sequence);
            std::swap(_seq_type, other._seq_type);

            return *this;
        }


        /**
         * \brief Compute and return the full assembling probability of this sequence (i.e., with all gene segments alignments).
         *
         * \param v_index Which aligned Variable gene segment to choose.
         * \param d_index Which aligned Diversity gene segment to choose.
         * \param j_index Which aligned Joining gene segment to choose.
         * \param action
         *
         * \return Full assembling probability. Returns 0 if wrong recombination function is used.
         */
        ///@{
        prob_t fullProbability(event_ind_t v_index, event_ind_t j_index) const {
            if (_recomb == VJ_RECOMB) {
                // P(Vi, Ji) * P(#dels | Vi) * P(V-J insertion seq) * P(#dels | Ji)
                return (matrix(0, 0)(v_index, j_index) *   // P(Vi & Ji)
                        matrix(1, v_index) *               // P(#dels | Vi)
                        matrix(2, 0) *                     // P(V-J insertion seq)
                        matrix(3, j_index))(0, 0);         // P(#dels | Ji)

            } else {
                return 0;
            }
        }

        prob_t fullProbability(event_ind_t v_index, event_ind_t d_index, event_ind_t j_index) const {
            if (_recomb == VDJ_RECOMB) {
                // P(Vi) * P(#dels | Vi) * P(V-D3' insertion seq) * P(D5'-D3' deletions | Di) * P(D5'-J insertion seq) * P(#dels | Ji) * P(Ji & Di)
                return (matrix(0, v_index) *      // P(Vi)
                        matrix(1, v_index) *      // P(#dels | Vi)
                        matrix(2, 0) *            // P(V-D3' insertion seq)
                        matrix(3, d_index) *      // P(D5'-D3' deletions | Di)
                        matrix(4, 0) *            // P(D5'-J insertion seq)
                        matrix(5, j_index) *      // P(#dels | Ji)
                        matrix(6, 0)(j_index, d_index))(0, 0);  // P(Ji & Di)
            } else {
                return 0;
            }
        }

        prob_t fullProbability(MAAGComputeProbAction action = SUM_PROBABILITY) const {
            // choose the max full probability from all possible recombinations of V(D)J gene segment indices
            if (_recomb == UNDEF_RECOMB) {
                return 0;
            }

            if (action == MAX_PROBABILITY) {
                prob_t max_prob = 0, cur_prob = 0;
                if (_recomb == VJ_RECOMB) {
                    for (event_ind_t v_index = 0; v_index < this->nVar(); ++v_index) {
                        for (event_ind_t j_index = 0; j_index < this->nJoi(); ++j_index) {
                            cur_prob = this->fullProbability(v_index, j_index);
                            if (cur_prob > max_prob) { max_prob = cur_prob; }
                        }
                    }
                } else {
                    for (event_ind_t v_index = 0; v_index < this->nVar(); ++v_index) {
                        for (event_ind_t d_index = 0; d_index < this->nDiv(); ++d_index) {
                            for (event_ind_t j_index = 0; j_index < this->nJoi(); ++j_index) {
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
                if (_recomb == VJ_RECOMB) {
                    for (event_ind_t v_index = 0; v_index < this->nVar(); ++v_index) {
                        for (event_ind_t j_index = 0; j_index < this->nJoi(); ++j_index) {
                            sum_prob += this->fullProbability(v_index, j_index);
                        }
                    }
                } else {
                    for (event_ind_t v_index = 0; v_index < this->nVar(); ++v_index) {
                        for (event_ind_t d_index = 0; d_index < this->nDiv(); ++d_index) {
                            for (event_ind_t j_index = 0; j_index < this->nJoi(); ++j_index) {
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
        event_ind_t nVar() const { return (_chain.size() == VJ_CHAIN_SIZE) ? this->nodeRows(VJ_VAR_JOI_GEN_I) : this->nodeSize(VDJ_VAR_GEN_I); }

        event_ind_t nJoi() const { return (_chain.size() == VJ_CHAIN_SIZE) ? this->nodeColumns(VJ_VAR_JOI_GEN_I) : this->nodeSize(VDJ_JOI_DEL_I); }

        event_ind_t nDiv() const { return (_chain.size() == VJ_CHAIN_SIZE) ? 0 : _chain[VDJ_DIV_DEL_I].size(); }
        ///@}


        /**
         * \brief Access to event indices, event probabilities and mismatches in the underlying matrices of the MAAG.
         *
         * \param node_i Index of the node.
         * \param mat_i Index of the matrix in the node.
         * \param row Which row to choose.
         * \param col Which column to choose.
         *
         * \return Event probability, event index or is mismatch. In the second case zero will be returned if no event chain matrix
         * is stored in this MAAG.
         */
        ///@{
        prob_t event_probability(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
            return (*this)(node_i, mat_i, row, col);
        }

        event_ind_t event_index(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
            return _events ? (*_events)(node_i, mat_i, row, col) : 0;
        };

        error_num_t errors(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
            return _errors ? (*_errors)(node_i, mat_i, row, col) : 0;
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
        bool save(const std::string &filepath) const {
            std::ofstream ofs(filepath);
            return this->save(ofs);
        }

        bool save(std::ostream &stream) const {
            return false;
        }

        bool load(const std::string &filepath) {
            std::ifstream ifs(filepath);
            return this->load(ifs);
        }

        bool load(std::istream &stream) {
            return false;
        }
        ///@}


        /**
         *
         */
        ///@{
        Recombination recombination() const { return _recomb; }

        bool is_vj() const { return _recomb == VJ_RECOMB; }

        bool is_vdj() const { return _recomb == VDJ_RECOMB; }

//        bool is_vd2j() const { return _recomb == VD2J_RECOMB; }
        ///@}


        bool has_events() const { return (bool) _events; }


        bool has_errors() const { return (bool) _errors; }

        prob_t error_prob() const { return _err_prob; }


        const std::string& sequence() const {
#ifndef DNDEBUG
            if (!_sequence) { throw(std::runtime_error("Access to a MAAG sequence when it's a nullptr!")); }
#endif
            return *_sequence;
        }


        seq_len_t n_poses() const { return _n_poses; }


        SequenceType sequence_type() const { return _seq_type; }

        dim_t rows(node_ind_t node_i) const { return this->nodeRows(node_i); }

        dim_t cols(node_ind_t node_i) const { return this->nodeColumns(node_i); }

    protected:

        pEventIndMMC _events;  /** Matrix of indices of events for each edge. */

        unique_ptr<seq_len_t[]> _seq_poses;  /** Vector of the initial clonotype sequence's positions for each vertex. */
        seq_len_t _n_poses;

        unique_ptr<sequence_t> _sequence;  /** Nucleotide or amino acid CDR3 sequence. */
        SequenceType _seq_type;

        pErrMMC _errors;  /** Matrix of number of errors for each scenario event position. */
        prob_t _err_prob;  /** Probability of a error. */

        Recombination _recomb;

//        MAAGMetadataPtr _metadata;

    };

}


//typedef MAAG<VJ_RECOMB, NO_METADATA, NO_ERRORS> MAAG_VJ_NM_NE;
//
//typedef MAAG<VJ_RECOMB, NO_METADATA, COMPUTE_ERRORS> MAAG_VJ_NM_ER;
//
//typedef MAAG<VJ_RECOMB, SAVE_METADATA, NO_ERRORS> MAAG_VJ_MD_NE;
//
//typedef MAAG<VJ_RECOMB, SAVE_METADATA, COMPUTE_ERRORS> MAAG_VJ_MD_ER;
//
//typedef MAAG<VDJ_RECOMB, NO_METADATA, NO_ERRORS> MAAG_VDJ_NM_NE;
//
//typedef MAAG<VDJ_RECOMB, NO_METADATA, COMPUTE_ERRORS> MAAG_VDJ_NM_ER;
//
//typedef MAAG<VDJ_RECOMB, SAVE_METADATA, NO_ERRORS> MAAG_VDJ_MD_NE;
//
//typedef MAAG<VDJ_RECOMB, SAVE_METADATA, COMPUTE_ERRORS> MAAG_VDJ_MD_ER;

#endif
