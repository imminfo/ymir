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
                  _insertions(other._insertions->clone()),
                  _insertions_rev(other._insertions_rev->clone()),
                  _max_ins_len(other._max_ins_len),
                  _max_ins_len_rev(other._max_ins_len_rev),
                  _ins_start(other._ins_start),
                  _ins_start_rev(other._ins_start_rev)
        {

        }


        MAAGaa(MAAGaa &&other)
                : MAAGBase(other)
        {
            _codons.swap(other._codons);
            _insertions.swap(other._insertions);
            _insertions_rev.swap(other._insertions_rev);
            std::swap(_max_ins_len, other._max_ins_len);
            std::swap(_max_ins_len_rev, other._max_ins_len_rev);
            std::swap(_ins_start, other._ins_start);
            std::swap(_ins_start_rev, other._ins_start_rev);
        }


        /**
         *
         */
        virtual ~MAAGaa() { }


        MAAGaa& operator= (const MAAGaa &other) {
            MAAGBase::operator=(other);
            _codons = other._codons;
            _insertions.reset(other._insertions->clone());
            _insertions_rev.reset(other._insertions_rev->clone());
            _max_ins_len = other._max_ins_len;
            _max_ins_len_rev = other._max_ins_len_rev;
            _ins_start = other._ins_start;
            _ins_start_rev = other._ins_start_rev;
            return *this;
        }


        MAAGaa& operator=(MAAGaa &&other) {
            MAAGBase::operator=(other);
            _codons.swap(other._codons);
            _insertions.swap(other._insertions);
            _insertions_rev.swap(other._insertions_rev);
            std::swap(_max_ins_len, other._max_ins_len);
            std::swap(_max_ins_len_rev, other._max_ins_len_rev);
            std::swap(_ins_start, other._ins_start);
            std::swap(_ins_start_rev, other._ins_start_rev);
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
#ifndef DNDEBUG
            if (_recomb != VJ_RECOMB) { return 0; }
#endif
            dim_t max_dim = std::max(this->nodeRows(2), this->nodeColumns(2));

            std::vector<prob_t> arr_prob1(max_dim);
            std::vector<prob_t> arr_prob2(max_dim);

            arr_prob1[0] = this->at(0, 0, v_index, j_index);
            for (dim_t i = 0; i < this->nodeColumns(1); ++i) {
                arr_prob2[i] = arr_prob1[0] * this->at(1, v_index, 0, i);
            }

            std::fill(arr_prob1.begin(), arr_prob1.end(), 0);

            for (dim_t row_i = 0; row_i < this->nodeRows(2); ++row_i) {
                for (dim_t col_i = 0; col_i < this->nodeColumns(2); ++col_i) {
                    int insertion_len = _seq_poses[this->nodeColumns(1) + col_i] - _seq_poses[row_i] - 1;
                    if (insertion_len >= 0 && insertion_len <= _max_ins_len) {
                        codon_hash left_codon, right_codon;
                        ((_seq_poses[row_i] / 3) == ((_seq_poses[row_i] + 1) / 3)) ? (_codons(0, v_index, 0, row_i))
                                                                                   : 63;

                        if (_seq_poses[row_i] == 0) {
                            left_codon = 63;
                        } else {
                            if (((_seq_poses[row_i] - 1) / 3) == ((_seq_poses[row_i]) / 3)) {
                                left_codon = _codons(0, v_index, 0, row_i);
                            } else {
                                left_codon = 63;
                            }
                        }

                        if (_seq_poses[this->nodeColumns(1) + col_i] == _sequence->size() * 3 + 1) {
                            right_codon = 63;
                        } else {
                            if (_seq_poses[this->nodeColumns(1) + col_i] <= 2) {
                                right_codon = 63;
                            } else {
                                if (((_seq_poses[this->nodeColumns(1) + col_i] - 2) / 3) ==
                                    ((_seq_poses[this->nodeColumns(1) + col_i] - 1) / 3)) {
                                    right_codon = _codons(1, j_index, col_i, 0);
                                } else {
                                    right_codon = 63;
                                }
                            }
                        }

                        prob_t prob_val = *(_ins_start + insertion_len)
                                          * _insertions->aaProbability(*_sequence,
                                                                       _seq_poses[row_i] + 1,
                                                                       _seq_poses[this->nodeColumns(1) + col_i] - 1,
                                                                       left_codon,
                                                                       right_codon);

                        arr_prob1[col_i] += prob_val * arr_prob2[row_i];
                    }
                }
            }

            arr_prob2[0] = 0;
            for (dim_t i = 0; i < this->nodeRows(3); ++i) {
                arr_prob2[0] += arr_prob1[i] * this->at(3, j_index, i, 0);
            }

            return arr_prob2[0];
        }

        prob_t fullProbability(event_ind_t v_index, event_ind_t d_index, event_ind_t j_index) {
#ifndef DNDEBUG
            if (_recomb != VDJ_RECOMB) { return 0; }
#endif
//                return (matrix(0, v_index) *      // P(Vi)
//                        matrix(1, v_index) *      // P(#dels | Vi)
//                        matrix(2, 0) *            // P(V-D3' insertion seq)
//                        matrix(3, d_index) *      // P(D5'-D3' deletions | Di)
//                        matrix(4, 0) *            // P(D5'-J insertion seq)
//                        matrix(5, j_index) *      // P(#dels | Ji)
//                        matrix(6, 0)(j_index, d_index))(0, 0);  // P(Ji & Di)

            dim_t max_dim = std::max(this->nodeRows(2), this->nodeColumns(2));
            max_dim = std::max(max_dim, this->nodeRows(4));
            max_dim = std::max(max_dim, this->nodeColumns(4));

            std::vector<prob_t> arr_prob1(max_dim);
            std::vector<prob_t> arr_prob2(max_dim);

            arr_prob1[0] = this->at(0, v_index, 0, 0);
            for (dim_t i = 0; i < this->nodeColumns(1); ++i) {
                arr_prob2[i] = arr_prob1[0] * this->at(1, v_index, 0, i);
            }

            std::fill(arr_prob1.begin(), arr_prob1.end(), 0);

            for (dim_t row_i = 0; row_i < this->nodeRows(2); ++row_i) {
                for (dim_t col_i = 0; col_i < this->nodeColumns(2); ++col_i) {
                    int insertion_len = _seq_poses[this->nodeColumns(1) + col_i] - _seq_poses[row_i] - 1;
                    if (insertion_len >= 0 && insertion_len <= _max_ins_len) {
                        codon_hash left_codon, right_codon;
                        ((_seq_poses[row_i] / 3) == ((_seq_poses[row_i] + 1) / 3)) ? (_codons(0, v_index, 0, row_i))
                                                                                   : 63;

                        if (_seq_poses[row_i] == 0) {
                            left_codon = 63;
                        } else {
                            if (((_seq_poses[row_i] - 1) / 3) == ((_seq_poses[row_i]) / 3)) {
                                left_codon = _codons(0, v_index, 0, row_i);
                            } else {
                                left_codon = 63;
                            }
                        }

                        if (_seq_poses[this->nodeColumns(1) + col_i] == _sequence->size() * 3 + 1) {
                            right_codon = 63;
                        } else {
                            if (_seq_poses[this->nodeColumns(1) + col_i] <= 2) {
                                right_codon = 63;
                            } else {
                                if (((_seq_poses[this->nodeColumns(1) + col_i] - 2) / 3) ==
                                    ((_seq_poses[this->nodeColumns(1) + col_i] - 1) / 3)) {
                                    right_codon = _codons(1, d_index, col_i, 0);
                                } else {
                                    right_codon = 63;
                                }
                            }
                        }

                        prob_t prob_val = *(_ins_start + insertion_len)
                                          * _insertions->aaProbability(*_sequence,
                                                                       _seq_poses[row_i] + 1,
                                                                       _seq_poses[this->nodeColumns(1) + col_i] - 1,
                                                                       left_codon,
                                                                       right_codon);

                        arr_prob1[col_i] += prob_val * arr_prob2[row_i];
                    }
                }
            }

            arr_prob2[0] = 0;
            for (dim_t i = 0; i < this->nodeRows(3); ++i) {
                arr_prob2[0] += arr_prob1[i] * this->at(3, j_index, i, 0);
            }

            return arr_prob2[0];
        }

        prob_t fullProbability(MAAGComputeProbAction action = SUM_PROBABILITY) {
            if (_recomb == UNDEF_RECOMB) {
                return 0;
            }

            // choose the max full probability from all possible recombinations of V(D)J gene segment indices
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


        codon_hash codon(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
            return _codons(node_i, mat_i, row, col);
        }


        void print() const {
            std::cout << *_sequence << std::endl;
            std::cout << matrix(0, 0).print() << std::endl;
            std::cout << matrix(1, 0).print() << std::endl;
            std::cout << matrix(2, 0).print() << std::endl;
            std::cout << matrix(3, 0).print() << std::endl;
            std::cout << "---------------------" << std::endl;
        }


        void printCodons() const {
            std::cout << "---------------------" << std::endl;
            std::cout << _codons.matrix(0, 0).print() << std::endl;
            std::cout << _codons.matrix(1, 0).print() << std::endl;
            std::cout << _codons.matrix(0, 1).print() << std::endl;
            std::cout << _codons.matrix(1, 1).print() << std::endl;
            std::cout << _codons.matrix(0, 2).print() << std::endl;
            std::cout << _codons.matrix(1, 2).print() << std::endl;
        }


    protected:

        CodonMMC _codons;
        unique_ptr<AbstractInsertionModel> _insertions, _insertions_rev;
        std::vector<prob_t>::const_iterator _ins_start, _ins_start_rev;
        seq_len_t _max_ins_len, _max_ins_len_rev;

    };

}


#endif //YMIR_MAAG_AA_H
