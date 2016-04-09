//
// Created by Vadim N. on 09/04/2016.
//

#ifndef YMIR_MAAG_NUC_H
#define YMIR_MAAG_NUC_H


#include "maag_base.h"


namespace ymir {

    class MAAGnuc;


    /**
     * \class MAAGnuc
     */
    class MAAGnuc : public MAAGBase {
    public:

        /**
         * \brief Default constructor.
         */
        MAAGnuc()
            : MAAGBase(NUCLEOTIDE),
              _events(nullptr),
              _errors(nullptr),
              _err_prob(0)
        {
            _values.reserve(1);
        }

        /**
         *
         */
        MAAGnuc(const MAAGnuc &other)
            : MAAGBase(other),
              _events(other._events ? new EventIndMMC(*other._events) : nullptr),
              _errors(other._errors ? new ErrMMC(*other._errors) : nullptr),
              _err_prob(other._err_prob)
        {
        }


        MAAGnuc(MAAGnuc &&other) {
            _events.swap(other._events);

            _errors.swap(other._errors);
            std::swap(_err_prob, other._err_prob);
        }


        /**
         *
         */
        virtual ~MAAGnuc() { }


        MAAGnuc& operator= (const MAAGnuc &other) {
            MAAGBase::operator=(other);

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

            return *this;
        }


        MAAGnuc& operator=(MAAGnuc &&other) {
            MAAGBase::operator=(other);

            _events.swap(other._events);

            _errors.swap(other._errors);
            std::swap(_err_prob, other._err_prob);

            return *this;
        }


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
        virtual prob_t fullProbability(event_ind_t v_index, event_ind_t j_index) const {
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

        virtual prob_t fullProbability(event_ind_t v_index, event_ind_t d_index, event_ind_t j_index) const {
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

        virtual prob_t fullProbability(MAAGComputeProbAction action = SUM_PROBABILITY) const {
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


    protected:

        prob_t _err_prob;  /** Probability of a error. */

        pEventIndMMC _events;  /** Matrix of indices of events for each edge. */
        pErrMMC _errors;  /** Matrix of number of errors for each scenario event position. */

    };

}


#endif //YMIR_MAAG_NUC_H
