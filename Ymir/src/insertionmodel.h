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

    class AbstractInsertionModel;
    class MonoNucInsertionModel;
    class DiNucInsertionModel;


    /**
    * \class AbstractInsertionModel
    *
    * \brief Class for representing VJ / VD / DJ insertions models - either a mono-nucleotide or a di-nucleotide model.
    */
    class AbstractInsertionModel {

    public:


        typedef unique_ptr<prob_t[]> prob_array_t;


        AbstractInsertionModel()
        {
        }


        AbstractInsertionModel(uint8_t size, prob_t err_prob = 0)
            : _size(size),
              _err_prob(err_prob),
              _arr(new prob_t[size])
        {
        }


        AbstractInsertionModel(uint8_t size, std::vector<prob_t>::const_iterator start, prob_t err_prob = 0)
            : AbstractInsertionModel(size, err_prob)
        {
            this->updateProbabilities(start);
        }


        AbstractInsertionModel(const AbstractInsertionModel &other)
                : _size(other._size),
                  _arr(new prob_t[other._size]),
                  _err_prob(other._err_prob)
        {
            this->updateProbabilities(other._arr.get());
        }


        AbstractInsertionModel& operator=(const AbstractInsertionModel &other)
        {
            _err_prob = other._err_prob;
            _size = other._size;
            if (_arr.get() != other._arr.get()) {
                _arr.reset(new prob_t[_size]);
                this->updateProbabilities(other._arr.get());
            }
            return *this;
        }


        virtual ~AbstractInsertionModel()
        {
        }


        /**
        * \brief Probability of the given nucleotide sequence.
        */
        ///@{
        virtual prob_t nucProbability(const std::string& sequence,
                                      char first_char = NULL_CHAR,
                                      bool with_errors = false) const = 0;

        virtual prob_t nucProbability(std::string::const_iterator start,
                                      seq_len_t sequence_len,
                                      char first_char = NULL_CHAR,
                                      bool with_errors = false) const = 0;

        virtual prob_t nucProbability(std::string::const_reverse_iterator start,
                                      seq_len_t sequence_len,
                                      char first_char = NULL_CHAR,
                                      bool with_errors = false) const = 0;
        ///@}


        /**
         *
         */
        ///@{
//        virtual prob_t aaProbability(const std::string& sequence,
//                                     seq_len_t first_nuc_pos,  // or an index of the codon? uint8_t first_aa_codons, second_aa_codons
//                                     seq_len_t last_nuc_pos,
//                                     char first_char = NULL_CHAR) const = 0;
//
//        virtual prob_t aaProbability(std::string::const_iterator start,
//                                     seq_len_t first_nuc_pos,
//                                     seq_len_t last_nuc_pos,
//                                     char first_char = NULL_CHAR) const = 0;
//
//        virtual prob_t aaProbability(std::string::const_reverse_iterator start,
//                                     seq_len_t first_nuc_pos,
//                                     seq_len_t last_nuc_pos,
//                                     char first_char = NULL_CHAR) const = 0;

        virtual prob_t aaProbability(const std::string& sequence, char first_char = NULL_CHAR) const
        {
            return this->aaProbability(sequence.cbegin(), sequence.size(), first_char);
        }

        virtual prob_t aaProbability(std::string::const_iterator start, seq_len_t sequence_len,
                                     char first_char = NULL_CHAR) const
        {
            throw(std::runtime_error("not implemented yet"));
        }

        virtual prob_t aaProbability(std::string::const_reverse_iterator start, seq_len_t sequence_len,
                                     char first_char = NULL_CHAR) const
        {
            throw(std::runtime_error("not implemented yet"));
        }
        ///@}


        /**
         *
         */
        virtual sequence_t generate(seq_len_t len, std::default_random_engine &rg, char first_char = NULL_CHAR, bool reverse = false) const = 0;


        /**
        * \name Event probability access.
        */
        ///@{
        prob_t operator[](uint8_t index) const { return _arr[index]; }
        ///@}


        prob_t err_prob() const { return _err_prob; }

    protected:

//        InsertionModelType _type;
        prob_array_t _arr;  /** Representation for a std::vector (#elements = 4) or a matrix (#elements = 16) with transition probabilities; rows and cols are for A-C-G-T (sequentially). */
        prob_t _err_prob;
        uint8_t _size;


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
        void updateProbabilities(std::vector<prob_t>::const_iterator start) {
            for (uint8_t i = 0; i < _size; ++i, ++start) {
                _arr[i] = *start;
            }
        }

        void updateProbabilities(prob_t *start) {
            for (uint8_t i = 0; i < _size; ++i) {
                _arr[i] = *(start + i);
            }
        }

        void updateProbabilities(std::vector<prob_t>::const_iterator start, prob_t err_prob) {
            this->updateProbabilities(start);
            _err_prob = err_prob;
        }

        void updateProbabilities(prob_t *start, prob_t err_prob) {
            this->updateProbabilities(start);
            _err_prob = err_prob;
        }
        ///@}

    };


    class MonoNucInsertionModel : public AbstractInsertionModel {
    public:

        MonoNucInsertionModel()
            : AbstractInsertionModel(4, 0)
        {
        }


        MonoNucInsertionModel(prob_t err_prob)
            : AbstractInsertionModel(4, err_prob)
        {
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
        prob_t aaProbability(const std::string& sequence,
                             char first_char = NULL_CHAR) const
        {
            return this->aaProbability(sequence.cbegin(), sequence.size(), first_char);
        }

        prob_t aaProbability(std::string::const_iterator start,
                             seq_len_t sequence_len,
                             char first_char = NULL_CHAR) const
        {
            throw(std::runtime_error("not implemented yet"));
        }

        prob_t aaProbability(std::string::const_reverse_iterator start,
                             seq_len_t sequence_len,
                             char first_char = NULL_CHAR) const
        {
            throw(std::runtime_error("not implemented yet"));
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

    protected:

    };


    class DiNucInsertionModel : public AbstractInsertionModel {
    public:


        DiNucInsertionModel()
            : AbstractInsertionModel(16, 0)
        {
        }


        DiNucInsertionModel(prob_t err_prob)
            : AbstractInsertionModel(16, err_prob)
        {
        }


        DiNucInsertionModel(std::vector<prob_t>::const_iterator start, prob_t err_prob = 0)
            : AbstractInsertionModel(16, err_prob)
        {
            this->updateProbabilities(start);
        }


        DiNucInsertionModel(const event_matrix_t& mat, prob_t err_prob = 0)
            : AbstractInsertionModel(16, err_prob)
        {
            this->updateProbabilitiesMatrix(mat);
        }


        DiNucInsertionModel(const DiNucInsertionModel &other)
            : AbstractInsertionModel(other)
        {
        }


        DiNucInsertionModel& operator=(const DiNucInsertionModel &other) {
            AbstractInsertionModel::operator=(other);
//            _err_prob = other._err_prob;
//            if (_arr.get() != other._arr.get()) {
//                _arr.reset(new prob_t[16]);
//                this->updateProbabilities(other._arr.get());
//            }
            return *this;
        }


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


        ///@{
        prob_t aaProbability(const std::string& sequence,
                             char first_char = NULL_CHAR) const
        {
            return this->aaProbability(sequence.cbegin(), sequence.size(), first_char);
        }

        prob_t aaProbability(std::string::const_iterator start,
                             seq_len_t sequence_len,
                             char first_char = NULL_CHAR) const
        {
            throw(std::runtime_error("not implemented yet"));
        }

        prob_t aaProbability(std::string::const_reverse_iterator start,
                             seq_len_t sequence_len,
                             char first_char = NULL_CHAR) const
        {
            throw(std::runtime_error("not implemented yet"));
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

    };


    // amino acid insertion model

}

#endif //YMIR_INSERTIONMODEL_H
