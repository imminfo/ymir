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
            codon_hash left_codon, right_codon;
            size_t left_pos, right_pos;
            int insertion_len;
            dim_t max_dim = std::max(this->rows(2), this->cols(2));

            std::vector<prob_t> arr_prob1(max_dim);
            std::vector<prob_t> arr_prob2(max_dim);

            // VJ genes + V deletions
            arr_prob1[0] = this->at(0, 0, v_index, j_index);
            for (dim_t i = 0; i < this->cols(1); ++i) {
                arr_prob2[i] = arr_prob1[0] * this->at(1, v_index, 0, i);
            }

            std::fill(arr_prob1.begin(), arr_prob1.end(), 0);

            // VJ insertions
            for (dim_t row_i = 0; row_i < this->rows(2); ++row_i) {
                for (dim_t col_i = 0; col_i < this->cols(2); ++col_i) {
                    left_pos = _seq_poses[row_i];
                    right_pos = _seq_poses[this->cols(1) + col_i];
                    insertion_len = right_pos - left_pos - 1;
                    if (insertion_len >= 0 && insertion_len <= _max_ins_len) {
                        if ((left_pos != 0) && ((left_pos - 1) / 3) == ((left_pos) / 3)) {
                            left_codon = _codons(0, v_index, 0, row_i);
                        } else {
                            left_codon = 63;
                        }

                        if ((right_pos != _sequence->size() * 3 + 1)
                            && (right_pos > 2)
                            && (((right_pos - 2) / 3) == ((right_pos - 1) / 3)))
                        {
                            right_codon = _codons(1, j_index, col_i, 0);
                        } else {
                            right_codon = 63;
                        }

                        prob_t prob_val = *(_ins_start + insertion_len)
                                          * _insertions->aaProbability(*_sequence,
                                                                       left_pos + 1,
                                                                       right_pos - 1,
                                                                       left_codon,
                                                                       right_codon);

                        arr_prob1[col_i] += prob_val * arr_prob2[row_i];
                    }
                }
            }

            arr_prob2[0] = 0;

            // J deletions
            for (dim_t i = 0; i < this->rows(3); ++i) {
                arr_prob2[0] += arr_prob1[i] * this->at(3, j_index, i, 0);
            }

            return arr_prob2[0];
        }

        prob_t fullProbability(event_ind_t v_index, event_ind_t d_index, event_ind_t j_index) {
#ifndef DNDEBUG
            if (_recomb != VDJ_RECOMB) { return 0; }
#endif

            codon_hash left_codon, right_codon, prev_codon;
            size_t left_pos, right_pos;
            int insertion_len;
            dim_t max_dim = std::max(this->rows(2), std::max(this->cols(2), this->cols(3)));

            std::vector<prob_t> arr_prob1(max_dim);
            std::vector<prob_t> arr_prob2(max_dim);

            // V gene + V deletions
            arr_prob1[0] = this->at(0, v_index, 0, 0);
            for (dim_t i = 0; i < this->cols(1); ++i) {
                arr_prob2[i] = arr_prob1[0] * this->at(1, v_index, 0, i);
            }

            std::fill(arr_prob1.begin(), arr_prob1.end(), 0);
            arr_prob1.resize(this->rows(3) * this->cols(3), 0);

            // VD insertions + D deletions
            for (dim_t row_i = 0; row_i < this->rows(2); ++row_i) {
                for (dim_t col_i = 0; col_i < this->cols(2); ++col_i) {
                    left_pos = _seq_poses[row_i];
                    right_pos = _seq_poses[this->cols(1) + col_i];
                    insertion_len = right_pos - left_pos - 1;
                    if (insertion_len >= 0 && insertion_len <= _max_ins_len) {
                        if ((left_pos != 0) && ((left_pos - 1) / 3) == (left_pos / 3)) {
                            left_codon = _codons(0, v_index, 0, row_i);

                            if (left_pos >= 1) {
                                prev_codon = _codons(0, v_index, 0, row_i - 1);
                            } else {
                                prev_codon = 63;
                            }
                        } else {
                            left_codon = 63;
                            prev_codon = 63;
                        }

                        for (dim_t d_col_i = 0; d_col_i < this->cols(3); ++d_col_i) {
                            if ((right_pos != _sequence->size() * 3 + 1)
                                && (right_pos > 2)
                                && (((right_pos - 2) / 3) == ((right_pos - 1) / 3)))
                            {
                                right_codon = _codons(1, d_index, col_i, d_col_i);
                            } else {
                                right_codon = 63;
                            }
//
//                            std::cout << "---" << std::endl;
//                            std::cout << (int)left_pos << std::endl;
//                            std::cout << (int)right_pos << std::endl;
//                            std::cout << (int)left_codon << std::endl;
//                            std::cout << (int)right_codon << std::endl;
//                            std::cout << (int)prev_codon << std::endl;

                            prob_t prob_val = *(_ins_start + insertion_len)
                                              * _insertions->aaProbability(*_sequence,
                                                                           left_pos + 1,
                                                                           right_pos - 1,
                                                                           left_codon,
                                                                           right_codon,
                                                                           prev_codon);

                            arr_prob1[col_i * this->cols(3) + d_col_i] += prob_val * arr_prob2[row_i] * this->at(3, d_index, col_i, d_col_i);
                        }
                    }
                }
            }

            std::fill(arr_prob2.begin(), arr_prob2.end(), 0);

            int start_dgen_nodes = this->cols(1) + this->cols(2),
                    start_dj_nodes = start_dgen_nodes + this->cols(3);

            // DJ insertions
            for (dim_t row_i = 0; row_i < this->rows(4); ++row_i) {
                for (dim_t col_i = 0; col_i < this->cols(4); ++col_i) {
                    left_pos = _seq_poses[start_dgen_nodes + row_i];
                    right_pos = _seq_poses[start_dj_nodes + col_i];
                    insertion_len = right_pos - left_pos - 1;
                    if (insertion_len >= 0 && insertion_len <= _max_ins_len_rev) {
                        if ((right_pos != _sequence->size() * 3 + 1)
                            && (right_pos > 2)
                            && (((right_pos - 2) / 3) == ((right_pos - 1) / 3)))
                        {
                            right_codon = _codons(3, j_index, col_i, 0);

                            if (right_pos < _sequence->size() * 3) {
                                prev_codon = _codons(3, j_index, col_i + 1, 0);
                            } else {
                                prev_codon = 63;
                            }
                        } else {
                            right_codon = 63;
                            prev_codon = 63;
                        }

                        for (dim_t d_row_i = 0; d_row_i < this->rows(3); ++d_row_i) {
                            if ((left_pos != 0) && ((left_pos - 1) / 3) == ((left_pos) / 3)) {
                                left_codon = _codons(2, d_index, d_row_i, row_i);
                            } else {
                                left_codon = 63;
                            }

                            prob_t prob_val = *(_ins_start_rev + insertion_len)
                                              * _insertions_rev->aaProbabilityRev(*_sequence,
                                                                                  right_pos - 1,
                                                                                  left_pos + 1,
                                                                                  right_codon,
                                                                                  left_codon,
                                                                                  prev_codon);

                            arr_prob2[col_i] += prob_val * arr_prob1[d_row_i * this->cols(3) + row_i];
                        }
                    }
                }
            }

            arr_prob1[0] = 0;

            // J deletions
            for (dim_t i = 0; i < this->rows(5); ++i) {
                arr_prob1[0] += arr_prob2[i] * this->at(5, j_index, i, 0);
            }

            // DJ genes
            return arr_prob1[0] * this->at(6, 0, j_index, d_index);
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
