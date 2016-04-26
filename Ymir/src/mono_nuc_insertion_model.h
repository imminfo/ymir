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
            this->make_aminoacid_probs();
        }


        MonoNucInsertionModel(const MonoNucInsertionModel &other)
            : AbstractInsertionModel(other),
              _aa_probs(other._aa_probs)
        {
        }


        MonoNucInsertionModel& operator=(const MonoNucInsertionModel &other) {
            AbstractInsertionModel::operator=(other);
//            _err_prob = other._err_prob;
//            if (_arr.get() != other._arr.get()) {
//                _arr.reset(new prob_t[4]);
//                this->updateProbabilities(other._arr.get());
//            }
            _aa_probs = other._aa_probs;
            return *this;
        }

        ~MonoNucInsertionModel()
        {
        }


        MonoNucInsertionModel* clone() const { return new MonoNucInsertionModel(*this); }


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
        virtual prob_t aaProbability(const sequence_t& sequence,
                                     seq_len_t first_nuc_pos,
                                     seq_len_t last_nuc_pos,
                                     codon_hash first_aa_codons,
                                     codon_hash last_aa_codons,
                                     char first_char = NULL_CHAR) const
        {
#ifndef DNDEBUG
            assert(first_nuc_pos <= 3*sequence.size() + 1);
            assert(last_nuc_pos <= 3*sequence.size() + 1);
#endif
            if (first_nuc_pos > last_nuc_pos) {
                return 1;
            }

            prob_t res = 0, nuc_prob;
            std::array<prob_t, 6> res_vec = {1, 1, 1, 1, 1, 1};

            prob_t arr_prob[5];
            arr_prob[0] = _arr[0];
            arr_prob[1] = _arr[1];
            arr_prob[2] = _arr[2];
            arr_prob[3] = _arr[3];
            arr_prob[4] = 0;

//            int aa_pos;
//            int codon_pos;
            std::array<int, 6> nuc_ids;

            if ((first_nuc_pos - 1) / 3 == (last_nuc_pos - 1) / 3) {
                std::bitset<6> bithash = first_aa_codons & last_aa_codons;

                for (int i = 0; i < 6; ++i) { res_vec[i] = bithash[5 - i]; }

                for (seq_len_t pos = first_nuc_pos; pos <= last_nuc_pos; ++pos) {
//                    aa_pos = (pos - 1) / 3;
//                    codon_pos = (pos - 1) % 3;
                    nuc_ids = CodonTable::table().which_nucl(sequence[(pos - 1) / 3], (pos - 1) % 3);
                    for (int i = 0; i < 6; ++i) {
                        res_vec[i] *= arr_prob[nuc_ids[i]];
                    }
                }

                for (int i = 0; i < 6; ++i) { res += res_vec[i]; }
            } else {
                prob_t tmp;

                // cycle through all amino acids
                std::bitset<6> bithash = first_aa_codons;
                for (int i = 0; i < 6; ++i) { res_vec[i] = bithash[5 - i]; }
                for (seq_len_t pos = first_nuc_pos; pos <= 3 * (1 + (first_nuc_pos - 1) / 3); ++pos) {
                    nuc_ids = CodonTable::table().which_nucl(sequence[(pos - 1) / 3], (pos - 1) % 3);
                    for (int i = 0; i < 6; ++i) {
                        res_vec[i] *= arr_prob[nuc_ids[i]];
                    }
                }
                for (int i = 0; i < 6; ++i) { res += res_vec[i]; }

                for (seq_len_t aa_pos = (1 + (first_nuc_pos - 1) / 3); aa_pos < ((last_nuc_pos - 1) / 3); ++aa_pos) {
                    res *= (*_aa_probs)[(sequence[aa_pos] << 8) + 63];

//                    res_vec = {1, 1, 1, 1, 1, 1};
//
//                    nuc_ids = CodonTable::table().which_nucl(sequence[aa_pos], 0);
//                    for (int i = 0; i < 6; ++i) { res_vec[i] *= arr_prob[nuc_ids[i]]; }
//
//                    nuc_ids = CodonTable::table().which_nucl(sequence[aa_pos], 1);
//                    for (int i = 0; i < 6; ++i) { res_vec[i] *= arr_prob[nuc_ids[i]]; }
//
//                    nuc_ids = CodonTable::table().which_nucl(sequence[aa_pos], 2);
//                    for (int i = 0; i < 6; ++i) { res_vec[i] *= arr_prob[nuc_ids[i]]; }
//
//                    tmp = 0;
//                    for (int i = 0; i < 6; ++i) { tmp += res_vec[i]; }
//                    res *= tmp;
                }

                bithash = last_aa_codons;
                for (int i = 0; i < 6; ++i) { res_vec[i] = bithash[5 - i]; }
                for (seq_len_t pos = 1 + 3 * ((last_nuc_pos - 1) / 3); pos <= last_nuc_pos; ++pos) {
                    nuc_ids = CodonTable::table().which_nucl(sequence[(pos - 1) / 3], (pos - 1) % 3);
                    for (int i = 0; i < 6; ++i) {
                        res_vec[i] *= arr_prob[nuc_ids[i]];
                    }
                }
                tmp = 0;
                for (int i = 0; i < 6; ++i) { tmp += res_vec[i]; }
                res *= tmp;
            }
            return res;
        }

        virtual prob_t aaProbabilityRev(const sequence_t &sequence,
                                        seq_len_t first_nuc_pos,
                                        seq_len_t last_nuc_pos,
                                        codon_hash first_aa_codons,
                                        codon_hash last_aa_codons,
                                        char first_char = NULL_CHAR) const
        {
            throw(std::runtime_error("aaProbabilityRev() is not implemented for MonoNucInsertionModel."));
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


    private:

        // aminoacid + codon hash -> prob
        // (aminoacid << 8) + codon_hash_value + start position
        // reverse or not???
        typedef shared_ptr<std::unordered_map<int16_t, prob_t>> shared_aa_ins_t;

        shared_aa_ins_t _aa_probs; // _aa_probs_frow, _aa_probs_back;
        // forward and reverse aminoacid probs


        void make_aminoacid_probs() {
            _aa_probs = std::make_shared<std::unordered_map<int16_t, prob_t>>();

            prob_t arr_prob[5];
            arr_prob[0] = _arr[0];
            arr_prob[1] = _arr[1];
            arr_prob[2] = _arr[2];
            arr_prob[3] = _arr[3];
            arr_prob[4] = 0;
            std::array<prob_t, 6> res_vec = {1, 1, 1, 1, 1, 1};

            for (auto it: CodonTable::table().aminoacids()) {
                if (it.first != '*') {
                    int16_t val = it.first << 8;
                    prob_t prob_val;
                    prob_t nuc_prob;

                    for (codon_hash hash_value = 0; hash_value <= 63; ++hash_value) {
//                    for (int start_pos = 0; start_pos <= 2; ++start_pos) {
                        // start position

                        std::bitset<6> bithash = hash_value;
                        for (int i = 0; i < 6; ++i) { res_vec[i] = bithash[5 - i]; }

                        auto nuc_ids = CodonTable::table().which_nucl(it.first, 0);
                        for (int i = 0; i < 6; ++i) { res_vec[i] *= arr_prob[nuc_ids[i]]; }

                        nuc_ids = CodonTable::table().which_nucl(it.first, 1);
                        for (int i = 0; i < 6; ++i) { res_vec[i] *= arr_prob[nuc_ids[i]]; }
//
                        nuc_ids = CodonTable::table().which_nucl(it.first, 2);
                        for (int i = 0; i < 6; ++i) { res_vec[i] *= arr_prob[nuc_ids[i]]; }

                        prob_t tmp = 0;
                        for (int i = 0; i < 6; ++i) { tmp += res_vec[i]; }

                        (*_aa_probs)[val + hash_value] = tmp;
//                        res *= tmp;
//                    }
                    }
                }
            }
        }

    };

}


#endif //YMIR_MONO_NUC_INSERTION_MODEL_H
