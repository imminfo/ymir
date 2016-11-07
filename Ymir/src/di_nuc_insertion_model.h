//
// Created by Vadim N. on 15/04/2016.
//

#ifndef YMIR_DI_NUC_INSERTION_MODEL_H
#define YMIR_DI_NUC_INSERTION_MODEL_H


#include "abstract_insertion_model.h"


namespace ymir {

    class DiNucInsertionModel : public AbstractInsertionModel {

        typedef std::bitset<6> bitset6;

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
         * \param prev_aa_codons
         */
        ///@{
        virtual prob_t aaProbability(const sequence_t& sequence,
                                     seq_len_t first_nuc_pos,
                                     seq_len_t last_nuc_pos,
                                     codon_hash first_aa_codons,
                                     codon_hash last_aa_codons,
                                     codon_hash prev_aa_codons = 1) const
        {

            // TODO: make _arr like that
//            prob_t arr_prob[16];
//            arr_prob[0] = _arr[0];
//            arr_prob[1] = _arr[1];
//            arr_prob[2] = _arr[2];
//            arr_prob[3] = _arr[3];
//            arr_prob[4] = 0;
//
//            arr_prob[5] = _arr[4];
//            arr_prob[6] = _arr[5];
//            arr_prob[7] = _arr[6];
//            arr_prob[8] = _arr[7];
//            arr_prob[9] = 0;
//
//            arr_prob[10] = _arr[8];
//            arr_prob[11] = _arr[9];
//            arr_prob[12] = _arr[10];
//            arr_prob[13] = _arr[11];
//            arr_prob[14] = 0;
//
//            arr_prob[15] = _arr[12];
//            arr_prob[16] = _arr[13];
//            arr_prob[17] = _arr[14];
//            arr_prob[18] = _arr[15];
//            arr_prob[19] = 0;
//
//            arr_prob[20] = 0;
//            arr_prob[21] = 0;
//            arr_prob[22] = 0;
//            arr_prob[23] = 0;
//            arr_prob[24] = 0;
//
//
//            std::array<prob_t, 6> res_vec;
//
//            if ((first_nuc_pos - 1) / 3 == (last_nuc_pos - 1) / 3) {
//                if (first_nuc_pos == 1) {
//                    res_vec = {.25, 0, 0, 0, 0, 0};
//                } else {
//                    (first_nuc_pos - 2) % 3
//                }
//
//                bitset6 bithash = first_aa_codons & last_aa_codons;
//
//                for (int i = 0; i < 6; ++i) { res_vec[i] *= bithash[5 - i]; }
//
//                for (seq_len_t pos = first_nuc_pos + 1; pos <= last_nuc_pos; ++pos) {
//                    nuc_ids = CodonTable::table().which_nucl(sequence[(pos - 1) / 3], (pos - 1) % 3);
//                    for (int i = 0; i < 6; ++i) {
//                        res_vec[i] *= arr_prob[nuc_ids[i]];
//                    }
//                }
//
//                for (int i = 0; i < 6; ++i) { res += res_vec[i]; }
//            } else {
//
//            }
        }

        virtual prob_t aaProbabilityRev(const sequence_t &sequence,
                                        seq_len_t first_nuc_pos,
                                        seq_len_t last_nuc_pos,
                                        codon_hash first_aa_codons,
                                        codon_hash last_aa_codons,
                                        codon_hash prev_aa_codons = 1) const
        {
            return 0;
        }
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

        typedef std::array<prob_t, 6> codon_prob_arr;

        // table of probs for each codon for each amino acid
        // (amino acid << 8) + (codon_hash_value << 2) + start (from) / end (to) position -> prob for each codon
        typedef shared_ptr<std::unordered_map<int16_t, codon_prob_arr>> shared_aa_ins_t;


        // table for probs of codons at the edge between two amino acids (i.e., P(x,y), where
        // x is a codon nuc on position 2 from the previous amino acid, and
        // y is a codon nuc on position 0 from the next amino acid on the chain)
        // for each amino acid for each hash
        // (amino acid, hash) x (amino acid, hash) -> prob for each codon
        // (amino acid1 << 8 + hash1) << 16 + (amino acid2 << 8 + hash2)
        typedef shared_ptr<std::unordered_map<int32_t, codon_prob_arr>> shared_aa_pair_t;

        // therefore the formula of probs for a pair of amino acids is the following:
        // R[0] = sum[i=0..5](shared_aa_ins_t(aa1+hash1+start_pos)[i] * shared_aa_pair_t(aa1+hash1+aa2+hash2)[i] * shared_aa_ins_t(aa2+hash2+to_pos)[0])
        // R[1] = sum[i=0..5](shared_aa_ins_t(aa1+hash1+start_pos)[i] * shared_aa_pair_t(aa1+hash1+aa2+hash2)[i] * shared_aa_ins_t(aa2+hash2+to_pos)[1])
        // ...


        shared_aa_ins_t _aa_probs_from;
        shared_aa_ins_t _aa_probs_to;
        shared_aa_pair_t _aa_probs_neis;


        /**
         *
         */
        ///@{
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
        ///@}


        /**
         *
         */
        void make_aminoacid_probs() {
            _aa_probs_forw = std::make_shared<std::unordered_map<char, codon_prob_arr>>();
            _aa_probs_back = std::make_shared<std::unordered_map<char, codon_prob_arr>>();
            _aa_probs_neis = std::make_shared<shared_aa_pair_t>();

            codon_prob_arr res_vec;

            for (auto it_prev: CodonTable::table().aminoacids()) {
                int16_t val_prev, val_next;
                auto codon_prev = CodonTable::table().codons(it_prev.first),
                        codon_next = CodonTable::table().codons(it_prev.first);

                for (auto it_next: CodonTable::table().aminoacids()) {
                    codon_next = CodonTable::table().codons(it_prev.next);

                    if (it_prev.first != '*' || it_next.first != "*") {
                        for (codon_hash hash_value_prev = 0; hash_value_prev <= 63; ++hash_value_prev) {
                            for (codon_hash hash_value_next = 0; hash_value_next <= 63; ++hash_value_next) {
                                val_prev = it_prev.first << 8;
                                val_next = it_next.first << 8;

                                res_vec.fill(0);

                                for (int i = 0; i < 6; ++i) { res_vec[i] = bithash[5 - i]; }

                                for (int i = 0; i < 6; ++i) {
                                    for (int j = 0; j < 6; ++j) {
                                        res_vec[i] *= (*this)(nuc_hash(codon_prev.codon()[2]), nuc_hash(codon_next.codon()[0]));
                                    }
                                }

                                (*_aa_probs_neis)[((val_prev + hash_value_prev) << 16) + (val_next + hash_value_next)] = res_vec;
                            }
                        }
                    }
                }
            }

            for (auto it: CodonTable::table().aminoacids()) {
                if (it.first != '*') {
                    int16_t val = it.first << 8;

                    res_vec.fill(0);

                    // there are two different modes to process the codons.
                    // is start pos == 0 then hash value should reflect the _previous_ codon's hash value.
                    // is start pos == 1 or 2 then hash values should reflect the _actual_ codon's hash value
                    // is start pos == 2 than all probs are 1s because _aa_probs_neis will reflect the
                    // transition probabilities to the next codon if needed
                    for (codon_hash hash_value = 0; hash_value <= 63; ++hash_value) {
                        bitset6 bithash = hash_value;

                        // start_pos == 0
                        for (int i = 0; i < 6; ++i) { res_vec[i] = bithash[5 - i]; }
                        for (int i = 0; i < 6; ++i) {
                            for (int j = 0; i < 6; ++i) {

                            }
                        }

                        (*_aa_probs_from)[val + (hash_value << 2) + 0] = res_vec;

                        // start_pos == 1 or 2
                        (*this)(nuc_hash(codon.codon()[0]), nuc_hash(codon.codon()[1]))
                        * (*this)(nuc_hash(codon.codon()[1]), nuc_hash(codon.codon()[2]));

                        (*_aa_probs_from)[val + (hash_value << 2) + 1] = res_vec;


                        (*_aa_probs_from)[val + (hash_value << 2) + 2] = res_vec;




                        for (int start_pos = 0; start_pos <= 2; ++start_pos) {
                            for (int i = 0; i < 6; ++i) { res_vec[i] = bithash[5 - i]; }

                            for (int pos = start_pos; pos <= 2; ++pos) {
                                auto nuc_ids = CodonTable::table().which_nucl(it.first, pos);
                                for (int i = 0; i < 6; ++i) { res_vec[i] *= arr_prob[nuc_ids[i]]; }
                            }

                            prob_t tmp = 0;
                            for (int i = 0; i < 6; ++i) { tmp += res_vec[i]; }
//                            std::cout << (int) start_pos << std::endl;
//                            std::cout << tmp << std::endl;

                            (*_aa_probs_from)[val + (hash_value << 2) + start_pos] = tmp;
                        }
                    }

                    for (codon_hash hash_value = 0; hash_value <= 63; ++hash_value) {
                        bitset6 bithash = hash_value;
                        for (int end_pos = 0; end_pos <= 2; ++end_pos) {
                            for (int i = 0; i < 6; ++i) { res_vec[i] = bithash[5 - i]; }

                            for (int pos = 0; pos <= end_pos; ++pos) {
                                auto nuc_ids = CodonTable::table().which_nucl(it.first, pos);
                                for (int i = 0; i < 6; ++i) { res_vec[i] *= arr_prob[nuc_ids[i]]; }
                            }

                            prob_t tmp = 0;
                            for (int i = 0; i < 6; ++i) { tmp += res_vec[i]; }

                            (*_aa_probs_to)[val + (hash_value << 2) + end_pos] = tmp;
                        }
                    }




                    auto codon = CodonTable::table().codons(it.first);

                    int i = 4;
                    codon_prob_arr prob_arr;
                    prob_arr.fill(0);
                    (*_aa_probs_forw)[it.first] = prob_arr;
                    (*_aa_probs_forw)[it.first][5] = (*this)(nuc_hash(codon.codon()[0]), nuc_hash(codon.codon()[1]))
                                                * (*this)(nuc_hash(codon.codon()[1]), nuc_hash(codon.codon()[2]));
                    while (codon.next()) {
                        (*_aa_probs_forw)[it.first][i] = (*this)(nuc_hash(codon.codon()[0]), nuc_hash(codon.codon()[1]))
                                                    * (*this)(nuc_hash(codon.codon()[1]), nuc_hash(codon.codon()[2]));
                        --i;
                    }

                    i = 4;
                    prob_arr.fill(0);
                    (*_aa_probs_back)[it.first] = prob_arr;
                    (*_aa_probs_back)[it.first][5] = (*this)(nuc_hash(codon.codon()[1]), nuc_hash(codon.codon()[0]))
                                                * (*this)(nuc_hash(codon.codon()[2]), nuc_hash(codon.codon()[1]));
                    while (codon.next()) {
                        (*_aa_probs_back)[it.first][i] = (*this)(nuc_hash(codon.codon()[1]), nuc_hash(codon.codon()[0]))
                                                    * (*this)(nuc_hash(codon.codon()[2]), nuc_hash(codon.codon()[1]));
                        --i;
                    }
                }
            }
        }

    };

}


#endif //YMIR_DI_NUC_INSERTION_MODEL_H
