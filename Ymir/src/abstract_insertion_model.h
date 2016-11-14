//
// Created by Vadim N. on 15/04/2016.
//

#ifndef YMIR_ABSTRACT_INSERTION_MODEL_H
#define YMIR_ABSTRACT_INSERTION_MODEL_H


#include <random>
#include <string>

#include "codon_alignment_vector.h"
#include "codon_table.h"
#include "types.h"
#include "tools.h"


namespace ymir {

    /**
    * \class AbstractInsertionModel
    *
    * \brief Class for representing VJ / VD / DJ insertions models - either a mono-nucleotide or a di-nucleotide model.
    */
    class AbstractInsertionModel {

    public:


        typedef std::bitset<6> bitset6;


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


        virtual AbstractInsertionModel* clone() const = 0;


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
        virtual prob_t aaProbability(const sequence_t& sequence,
                                     seq_len_t first_nuc_pos,
                                     seq_len_t last_nuc_pos,
                                     codon_hash first_aa_codons,
                                     codon_hash last_aa_codons,
                                     codon_hash prev_aa_codons = 1) const = 0;

        virtual prob_t aaProbabilityRev(const sequence_t &sequence,
                                        seq_len_t first_nuc_pos,
                                        seq_len_t last_nuc_pos,
                                        codon_hash first_aa_codons,
                                        codon_hash last_aa_codons,
                                        codon_hash prev_aa_codons = 1) const = 0;
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


        virtual void make_aminoacid_probs() = 0;

    };

}

#endif //YMIR_ABSTRACT_INSERTION_MODEL_H
