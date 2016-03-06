//
// Created by Vadim N. on 21/06/2015.
//

#ifndef YMIR_INSERTIONMODEL_H
#define YMIR_INSERTIONMODEL_H

#include <random>
#include <string>

#include "types.h"
#include "tools.h"


namespace ymir {

    /**
    * \class InsertionModel
    *
    * \brief Class for representing VJ / VD / DJ insertions models - either a mono-nucleotide or a di-nucleotide model.
    */
    class InsertionModel {

    public:


        typedef unique_ptr<prob_t[]> prob_array_t;


        InsertionModel(InsertionModelType mt, prob_t err_prob = 0)
            : _type(mt),
              _arr(new prob_t[mt == MONO_NUCLEOTIDE ? 4 : 16]),
              _err_prob(err_prob)
        {
        }


        InsertionModel(InsertionModelType mt, std::vector<prob_t>::const_iterator start, prob_t err_prob = 0)
            : InsertionModel(mt, err_prob)
        {
            this->updateProbabilities(start);
        }


        InsertionModel(const event_matrix_t& mat, prob_t err_prob = 0)
            : InsertionModel(DI_NUCLEOTIDE, err_prob)
        {
            this->updateProbabilities(mat);
        }


        InsertionModel(const InsertionModel &other) 
            : _type(other._type), 
              _arr(new prob_t[other._type == MONO_NUCLEOTIDE ? 4 : 16]),
              _err_prob(other._err_prob)
        {
            this->updateProbabilities(other._arr.get());
        }


        InsertionModel& operator=(const InsertionModel &other) {
            _err_prob = other._err_prob;
            if (_arr.get() != other._arr.get()) {
                _type = other._type;
                this->initProbabilities();
                this->updateProbabilities(other._arr.get());
            }
            return *this;
        }


        ~InsertionModel()
        {
        }


        void initProbabilities() {
            _arr.reset(new prob_t[_type == MONO_NUCLEOTIDE ? 4 : 16]);
            _err_prob = 0;
        }


        /**
        * \name Update chain probabilities.
        *
        * \brief Update internal std::vector of nucleotide joint probabilities.
        *
        * \param start Iterator to std::vector, in which each 4 elements are rows in nucleotide joint probability matrix, i.e.
        * row for prev-A-nucleotide, row for prev-C-nucleotide and so on.
        * \param mat Matrix with nucleotide joint probabilities with rows for prev nucleotide and columns
        * for next nucleotide.
        */
        ///@{
        void updateProbabilities(std::vector<prob_t>::const_iterator start, prob_t err_prob = 0) {
            uint8_t last = _type == MONO_NUCLEOTIDE ? 4 : 16;
            for (uint8_t i = 0; i < last; ++i, ++start) {
                _arr[i] = *start;
            }
            _err_prob = err_prob;
        }

        void updateProbabilities(prob_t *start, prob_t err_prob = 0) {
            uint8_t last = _type == MONO_NUCLEOTIDE ? 4 : 16;
            for (uint8_t i = 0; i < last; ++i) {
                _arr[i] = *(start + i);
            }
            _err_prob = err_prob;
        }

        void updateProbabilities(const event_matrix_t& mat, prob_t err_prob = 0) {
#ifndef DNDEBUG
            if (_type == MONO_NUCLEOTIDE) { throw(std::runtime_error("Can't initialise the mono-nucleotide insertion model with a matrix!")); }
#endif
            for (uint8_t i = 0; i < 4; ++i) {
                for (uint8_t j = 0; j < 4; ++j) {
                    _arr[4*i + j] = mat(i, j);
                }
            }

            _err_prob = err_prob;
        }
        ///@}


        /**
        * \brief Probability of the given nucleotide sequence.
        */
        ///@{
        prob_t nucProbability(const std::string& sequence, char first_char = NULL_CHAR) const {
            return this->nucProbability<std::string::const_iterator>(sequence.cbegin(), sequence.size(), first_char);
        }

        template <class Iterator>
        prob_t nucProbability(Iterator start, seq_len_t sequence_len, char first_char = NULL_CHAR) const {
            prob_t res = 1;

            if (sequence_len) {
                if (_type == MONO_NUCLEOTIDE) {
                    auto tmp = start;
                    for (seq_len_t i = 0; i < sequence_len; ++i, ++start) {
                        res *= _arr[nuc_hash(*start)];
                    }
//                    if (std::isnan(res)) {
//                        cout << "NAN res!!" << endl;
//                        res = 1;
//                        tmp = start;
//                        for (seq_len_t i = 0; i < sequence_len; ++i, ++tmp) {
//                            cout << "char:" << (*tmp) << " " << res << " -> ";
//                            res *= _arr[nuc_hash(*tmp)];
//                            cout << res << endl;
//                        }
//                    }
                } else {
                    auto tmp1 = start, tmp2 = start + 1;
                    auto next = start + 1;
                    res = (first_char == NULL_CHAR) ? .25 : (*this)(nuc_hash(first_char), nuc_hash(*start));
                    for (seq_len_t i = 1; i < sequence_len; ++i, ++start, ++next) {
                        res *= (*this)(nuc_hash(*start), nuc_hash(*next));
                    }
//                    if (std::isnan(res)) {
//                        tmp1 = start;
//                        tmp2 = start + 1;
//                        cout << "NAN res!!" << endl;
//                        res = (first_char == NULL_CHAR) ? .25 : (*this)(nuc_hash(first_char), nuc_hash(*tmp1));
//                        cout << "first char: " << first_char << endl;
//                        cout << "res1: " << res << endl;
//                        cout << "is null char: " << (first_char == NULL_CHAR) << endl;
//                        cout << "hash: " << ((*this)(nuc_hash(first_char), nuc_hash(*tmp1))) << endl;
//                        cout << "res2: " << res << endl;
//                        for (seq_len_t i = 1; i < sequence_len; ++i, ++tmp1, ++tmp2) {
//                            cout << "char:" << (*tmp1) << ":" << (*tmp2) << " " << res << " -> ";
//                            res *= (*this)(nuc_hash(*tmp1), nuc_hash(*tmp2));
//                            cout << res << endl;
//                        }
//                    }
                }
            }

            return res;
        }
        ///@}


        /**
         *
         */
        ///@{
        prob_t aaProbability(const std::string& sequence, char first_char = NULL_CHAR) const {
            return this->aaProbability(sequence.cbegin(), sequence.size(), first_char);
        }

        prob_t aaProbability(std::string::const_iterator start, seq_len_t sequence_len, char first_char = NULL_CHAR) const {
            return 0;
        }
        ///@}


        /**
         *
         */
        std::string generate(seq_len_t len, std::default_random_engine &rg, char first_char = NULL_CHAR, bool reverse = false) const {
            std::string res = "";
            if (len) {
                if (_type == MONO_NUCLEOTIDE) {
                    std::discrete_distribution<int> distr;
                    distr = std::discrete_distribution<int>(_arr.get(), _arr.get() + 4);

                    for (seq_len_t i = 0; i < len; ++i) {
                        res += inv_nuc_hash(distr(rg));
                    }
                } else {
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
            }

            return res;
        }


        /**
        * \name Event probability access.
        */
        ///@{
        prob_t operator[](uint8_t index) const { return _arr[index]; }
        prob_t operator()(uint8_t row, uint8_t col) const { return _arr[4*row + col]; }
        ///@}

    protected:

        InsertionModelType _type;
        prob_array_t _arr;  /** Representation for a std::vector (#elements = 4) or a matrix (#elements = 16) with transition probabilities; rows and cols are for A-C-G-T (sequentially). */
        prob_t _err_prob;


        InsertionModel() 
            : _arr(new prob_t[4])
        {
            _type = MONO_NUCLEOTIDE;
        }

    };
}

#endif //YMIR_INSERTIONMODEL_H
