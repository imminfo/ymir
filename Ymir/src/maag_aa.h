//
// Created by Vadim N. on 09/04/2016.
//

#ifndef YMIR_MAAG_AA_H
#define YMIR_MAAG_AA_H


#include "maag_base.h"


namespace ymir {

    class MAAGaa;


    /**
     * \class MAAGaa
     */
    class MAAGaa : public MAAGBase {

        friend class MAAGBuilder;
        friend class MAAGForwardBackwardAlgorithm;

    public:

        /**
         * \brief Default constructor.
         */
        MAAGaa()
                : MAAGBase(AMINOACID) //, codons, insertions
        {
            _values.reserve(1);
        }

        /**
         *
         */
        MAAGaa(const MAAGaa &other)
                : MAAGBase(other),
                  _codons(other._codons),
                  _insertions(other._insertions->clone())
        {

        }


        MAAGaa(MAAGaa &&other)
                : MAAGBase(other)
        {
            _codons.swap(other._codons);
            _insertions.swap(other._insertions);
        }


        /**
         *
         */
        virtual ~MAAGaa() { }


        MAAGaa& operator= (const MAAGaa &other) {
            MAAGBase::operator=(other);
            _codons = other._codons;
            _insertions.reset(other._insertions->clone());
            return *this;
        }


        MAAGaa& operator=(MAAGaa &&other) {
            MAAGBase::operator=(other);
            _codons.swap(other._codons);
            _insertions.swap(other._insertions);
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
        prob_t fullProbability(event_ind_t v_index, event_ind_t j_index) {
            if (_recomb == VJ_RECOMB) {
                // build insertion matrix
                this->fill(2, 0, 0);
                for (dim_t row_i = 0; row_i < this->nodeRows(2); ++row_i) {
                    for (dim_t col_i = 0; col_i < this->nodeColumns(2); ++col_i) {
                        if (_seq_poses[this->nodeColumns(1) + col_i] >= _seq_poses[row_i]) {
                            this->at(2, 0, row_i, col_i) = _insertions->aaProbability(*_sequence,
                                                                                      _seq_poses[row_i],
                                                                                      _seq_poses[this->nodeColumns(1) + col_i],
                                                                                      _codons(0, v_index, 0, row_i),
                                                                                      _codons(1, j_index, col_i, 0));
                        }
                    }
                }

                // compute the full probability
                // P(Vi, Ji) * P(#dels | Vi) * P(V-J insertion seq) * P(#dels | Ji)
                return (matrix(0, 0)(v_index, j_index) *   // P(Vi & Ji)
                        matrix(1, v_index) *               // P(#dels | Vi)
                        matrix(2, 0) *                     // P(V-J insertion seq)
                        matrix(3, j_index))(0, 0);         // P(#dels | Ji)

            } else {
                return 0;
            }
        }

        prob_t fullProbability(event_ind_t v_index, event_ind_t d_index, event_ind_t j_index) {
//            if (_recomb == VDJ_RECOMB) {
//                // P(Vi) * P(#dels | Vi) * P(V-D3' insertion seq) * P(D5'-D3' deletions | Di) * P(D5'-J insertion seq) * P(#dels | Ji) * P(Ji & Di)
//                return (matrix(0, v_index) *      // P(Vi)
//                        matrix(1, v_index) *      // P(#dels | Vi)
//                        matrix(2, 0) *            // P(V-D3' insertion seq)
//                        matrix(3, d_index) *      // P(D5'-D3' deletions | Di)
//                        matrix(4, 0) *            // P(D5'-J insertion seq)
//                        matrix(5, j_index) *      // P(#dels | Ji)
//                        matrix(6, 0)(j_index, d_index))(0, 0);  // P(Ji & Di)
//            } else {
//                return 0;
//            }
            return 0;
        }

        prob_t fullProbability(MAAGComputeProbAction action = SUM_PROBABILITY) {
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

        CodonMMC _codons;
        unique_ptr<AbstractInsertionModel> _insertions;

    };

}


#endif //YMIR_MAAG_AA_H
