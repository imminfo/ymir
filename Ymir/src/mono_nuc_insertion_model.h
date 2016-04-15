//
// Created by Vadim N. on 15/04/2016.
//

#ifndef YMIR_MONO_NUC_INSERTION_MODEL_H
#define YMIR_MONO_NUC_INSERTION_MODEL_H


#include "abstract_insertion_model.h"


namespace ymir {

    class MonoNucInsertionModel : public AbstractInsertionModel {

    public:

        MonoNucInsertionModel()
            : AbstractInsertionModel(4, 0)
        {
            this->make_aminoacid_probs();
        }


        MonoNucInsertionModel(prob_t err_prob)
            : AbstractInsertionModel(4, err_prob)
        {
            this->make_aminoacid_probs();
        }


        MonoNucInsertionModel(std::vector<prob_t>::const_iterator start, prob_t err_prob = 0)
            : AbstractInsertionModel(4, err_prob)
        {
            this->updateProbabilities(start);
        }


        MonoNucInsertionModel(const MonoNucInsertionModel &other)
            : AbstractInsertionModel(other)
        {
        }


        MonoNucInsertionModel& operator=(const MonoNucInsertionModel &other) {
            AbstractInsertionModel::operator=(other);
//            _err_prob = other._err_prob;
//            if (_arr.get() != other._arr.get()) {
//                _arr.reset(new prob_t[4]);
//                this->updateProbabilities(other._arr.get());
//            }
            return *this;
        }


        prob_t nucProbability(const std::string& sequence,
                              char first_char = NULL_CHAR,
                              bool with_errors = false) const
        {
            return this->nucProbability(sequence.cbegin(), sequence.size(), first_char, with_errors);
        }

        prob_t nucProbability(std::string::const_iterator start,
                              seq_len_t sequence_len,
                              char first_char = NULL_CHAR,
                              bool with_errors = false) const
        {
            prob_t res = 1;

            if (sequence_len) {
                auto tmp = start;
                if (!with_errors) {
                    for (seq_len_t i = 0; i < sequence_len; ++i, ++start) {
                        res *= _arr[nuc_hash(*start)];
                    }
                } else {
                    for (seq_len_t i = 0; i < sequence_len; ++i, ++start) {
                        res *= (_err_prob + (1 - _err_prob) * _arr[nuc_hash(*start)]);
                    }
                }
            }

            return res;
        }

        prob_t nucProbability(std::string::const_reverse_iterator start,
                              seq_len_t sequence_len,
                              char first_char = NULL_CHAR,
                              bool with_errors = false) const
        {
            prob_t res = 1;

            if (sequence_len) {
                auto tmp = start;
                if (!with_errors) {
                    for (seq_len_t i = 0; i < sequence_len; ++i, ++start) {
                        res *= _arr[nuc_hash(*start)];
                    }
                } else {
                    for (seq_len_t i = 0; i < sequence_len; ++i, ++start) {
                        res *= (_err_prob + (1 - _err_prob) * _arr[nuc_hash(*start)]);
                    }
                }
            }

            return res;
        }


        /**
         *
         */
        ///@{
        virtual prob_t aaProbability(const std::string& sequence,
                                     seq_len_t first_nuc_pos,
                                     seq_len_t last_nuc_pos,
                                     codon_hash first_aa_codons,
                                     codon_hash last_aa_codons,
                                     char first_char = NULL_CHAR) const
        {

        }

        virtual prob_t aaProbability(std::string::const_iterator start,
                                     seq_len_t sequence_len,
                                     seq_len_t first_nuc_pos,
                                     seq_len_t last_nuc_pos,
                                     codon_hash first_aa_codons,
                                     codon_hash last_aa_codons,
                                     char first_char = NULL_CHAR) const
        {

        }

        virtual prob_t aaProbability(std::string::const_reverse_iterator start,
                                     seq_len_t sequence_len,
                                     seq_len_t first_nuc_pos,
                                     seq_len_t last_nuc_pos,
                                     codon_hash first_aa_codons,
                                     codon_hash last_aa_codons,
                                     char first_char = NULL_CHAR) const
        {

        }
        ///@}


        sequence_t generate(seq_len_t len, std::default_random_engine &rg, char first_char = NULL_CHAR, bool reverse = false) const {
            std::string res = "";
            if (len) {
                std::discrete_distribution<int> distr;
                distr = std::discrete_distribution<int>(_arr.get(), _arr.get() + 4);

                for (seq_len_t i = 0; i < len; ++i) {
                    res += inv_nuc_hash(distr(rg));
                }
            }

            return res;
        }


        void make_aminoacid_probs() {

        }

    };

}


#endif //YMIR_MONO_NUC_INSERTION_MODEL_H
