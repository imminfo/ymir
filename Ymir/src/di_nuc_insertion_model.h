//
// Created by Vadim N. on 15/04/2016.
//

#ifndef YMIR_DI_NUC_INSERTION_MODEL_H
#define YMIR_DI_NUC_INSERTION_MODEL_H


#include "abstract_insertion_model.h"


namespace ymir {

    class DiNucInsertionModel : public AbstractInsertionModel {
    public:


        DiNucInsertionModel()
            : AbstractInsertionModel(16, 0)
        {
            this->make_aminoacid_probs();
        }


        DiNucInsertionModel(prob_t err_prob)
            : AbstractInsertionModel(16, err_prob)
        {
            this->make_aminoacid_probs();
        }


        DiNucInsertionModel(std::vector<prob_t>::const_iterator start, prob_t err_prob = 0)
            : AbstractInsertionModel(16, err_prob)
        {
            this->updateProbabilities(start);
            this->make_aminoacid_probs();
        }


        DiNucInsertionModel(const event_matrix_t& mat, prob_t err_prob = 0)
            : AbstractInsertionModel(16, err_prob)
        {
            this->updateProbabilitiesMatrix(mat);
            this->make_aminoacid_probs();
        }


        DiNucInsertionModel(const DiNucInsertionModel &other)
            : AbstractInsertionModel(other),
              _aa_probs(other._aa_probs)
        {
        }


        DiNucInsertionModel& operator=(const DiNucInsertionModel &other) {
            AbstractInsertionModel::operator=(other);
//            _err_prob = other._err_prob;
//            if (_arr.get() != other._arr.get()) {
//                _arr.reset(new prob_t[16]);
//                this->updateProbabilities(other._arr.get());
//            }
            _aa_probs = other._aa_probs;
            return *this;
        }


        ~DiNucInsertionModel()
        {
        }


        virtual DiNucInsertionModel* clone() const { return new DiNucInsertionModel(*this); }


        ///@{
        prob_t nucProbability(const std::string& sequence, char first_char = NULL_CHAR, bool with_errors = false) const {
            return this->nucProbability(sequence.cbegin(), sequence.size(), first_char, with_errors);
        }

        prob_t nucProbability(std::string::const_iterator start, seq_len_t sequence_len, char first_char = NULL_CHAR, bool with_errors = false) const {
            prob_t res = 1;

            if (sequence_len) {
                auto tmp1 = start, tmp2 = start + 1;
                auto next = start + 1;
                if (!with_errors) {
                    res = (first_char == NULL_CHAR) ? .25 : (*this)(nuc_hash(first_char), nuc_hash(*start));
                    for (seq_len_t i = 1; i < sequence_len; ++i, ++start, ++next) {
                        res *= (*this)(nuc_hash(*start), nuc_hash(*next));
                    }
                } else {
                    prob_t err2 = _err_prob * _err_prob;
                    res = (first_char == NULL_CHAR) ? .25 : (err2 + (1 - err2) * (*this)(nuc_hash(first_char), nuc_hash(*start)));
                    for (seq_len_t i = 1; i < sequence_len; ++i, ++start, ++next) {
                        res *= (err2 + (1 - err2) * (*this)(nuc_hash(*start), nuc_hash(*next)));
                    }
                }
            }

            return res;
        }

        prob_t nucProbability(std::string::const_reverse_iterator start, seq_len_t sequence_len, char first_char = NULL_CHAR, bool with_errors = false) const {
            prob_t res = 1;

            if (sequence_len) {
                auto tmp1 = start, tmp2 = start + 1;
                auto next = start + 1;
                if (!with_errors) {
                    res = (first_char == NULL_CHAR) ? .25 : (*this)(nuc_hash(first_char), nuc_hash(*start));
                    for (seq_len_t i = 1; i < sequence_len; ++i, ++start, ++next) {
                        res *= (*this)(nuc_hash(*start), nuc_hash(*next));
                    }
                } else {
                    prob_t err2 = _err_prob * _err_prob;
                    res = (first_char == NULL_CHAR) ? .25 : (err2 + (1 - err2) * (*this)(nuc_hash(first_char), nuc_hash(*start)));
                    for (seq_len_t i = 1; i < sequence_len; ++i, ++start, ++next) {
                        res *= (err2 + (1 - err2) * (*this)(nuc_hash(*start), nuc_hash(*next)));
                    }
                }
            }

            return res;
        }
        ///@}


        /**
         *
         */
        ///@{
        virtual prob_t aaProbability(const sequence_t& sequence,
                                     seq_len_t first_nuc_pos,
                                     seq_len_t last_nuc_pos,
                                     codon_hash first_aa_codons,
                                     codon_hash last_aa_codons,
                                     char first_char = NULL_CHAR) const {}

        virtual prob_t aaProbabilityRev(const sequence_t &sequence,
                                        seq_len_t first_nuc_pos,
                                        seq_len_t last_nuc_pos,
                                        codon_hash first_aa_codons,
                                        codon_hash last_aa_codons,
                                        char first_char = NULL_CHAR) const {}
        ///@}


        sequence_t generate(seq_len_t len, std::default_random_engine &rg, char first_char = NULL_CHAR, bool reverse = false) const {
            std::string res = "";
            if (len) {
                std::discrete_distribution<int> distrs[] = {
                        std::discrete_distribution<int>(_arr.get(), _arr.get() + 4),
                        std::discrete_distribution<int>(_arr.get() + 4, _arr.get() + 8),
                        std::discrete_distribution<int>(_arr.get() + 8, _arr.get() + 12),
                        std::discrete_distribution<int>(_arr.get() + 12, _arr.get() + 16)
                };

                if (first_char == NULL_CHAR) {
                    res = inv_nuc_hash(distrs[std::discrete_distribution<int>{.25, .25, .25, .25}(rg)](rg));
                } else {
                    res = inv_nuc_hash(distrs[nuc_hash(first_char)](rg));
                }
                for (seq_len_t i = 1; i < len; ++i) {
                    res += inv_nuc_hash(distrs[nuc_hash(res[i - 1])](rg));
                }

                if (reverse) { std::reverse(res.begin(), res.end()); }
            }

            return res;
        }


        prob_t operator()(uint8_t row, uint8_t col) const { return _arr[4*row + col]; }

    protected:


        typedef shared_ptr<std::unordered_map<char, prob_t>> shared_aa_ins_t;

        shared_aa_ins_t _aa_probs;


        void updateProbabilitiesMatrix(const event_matrix_t& mat) {
            for (uint8_t i = 0; i < 4; ++i) {
                for (uint8_t j = 0; j < 4; ++j) {
                    _arr[4*i + j] = mat(i, j);
                }
            }
        }

        void updateProbabilitiesMatrix(const event_matrix_t& mat, prob_t err_prob) {
            this->updateProbabilitiesMatrix(mat);
            _err_prob = err_prob;
        }


        void make_aminoacid_probs() {

        }

    };

}


#endif //YMIR_DI_NUC_INSERTION_MODEL_H
