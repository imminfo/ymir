//
// Created by Vadim N. on 21/06/2015.
//

#ifndef YMIR_INSERTIONMODEL_H
#define YMIR_INSERTIONMODEL_H

#include "string"

#include "types.h"


namespace ymir {

    /**
    * \class InsertionModel
    *
    * \brief Class for representing VJ / VD / DJ insertions models - either a mono-nucleotide or a di-nucleotide model.
    */
    class InsertionModel {
    public:

        InsertionModel(INSERTION_MODEL_TYPE mt) {
            this->_type = mt;
            this->initProbabilities();
        }


        InsertionModel(INSERTION_MODEL_TYPE mt, vector<prob_t>::const_iterator start) : InsertionModel(mt) {
            this->updateProbabilities(start);
        }


        InsertionModel(const event_matrix_t& mat) : InsertionModel(DiNucleotide) {
            this->updateProbabilities(mat);
        }


        ~InsertionModel() {
            delete [] this->_arr;
        }


        void initProbabilities() {
            if (_type == MonoNucleotide) {
                this->_arr = new prob_t[4];
            } else {
                this->_arr = new prob_t[16];
            }
        }


        /**
        * \brief Probability of the given nucleotide sequence.
        */
        prob_t nucProbability(const string& sequence, char first_char = NULL_CHAR) const {
            return this->nucProbability(sequence.cbegin(), sequence.size(), first_char);
        }

        prob_t nucProbability(string::const_iterator start, seq_len_t sequence_len, char first_char = NULL_CHAR) const {
            prob_t res = 1;

            if (sequence_len) {
                if (_type == MonoNucleotide) {
                    auto tmp = start;
                    for (seq_len_t i = 0; i < sequence_len; ++i, ++start) {
                        res *= _arr[nuc_hash(*start)];
                    }
                    if (isnan(res)) {
                        cout << "NAN res!!" << endl;
                        res = 1;
                        for (seq_len_t i = 0; i < sequence_len; ++i, ++tmp) {
                            cout << "char:" << (*tmp) << " " << res << " -> ";
                            res *= _arr[nuc_hash(*tmp)];
                            cout << res << endl;
                        }
                    }
                } else {
                    auto tmp1 = start, tmp2 = start + 1;
                    string::const_iterator next = start + 1;
                    res = (first_char == NULL_CHAR) ? .25 : (*this)(nuc_hash(first_char), nuc_hash(*start));
                    for (seq_len_t i = 1; i < sequence_len; ++i, ++start, ++next) {
                        res *= (*this)(nuc_hash(*start), nuc_hash(*next));
                    }
//                    if (isnan(res)) {
//                        cout << "NAN res!!" << endl;
//                        res = (first_char == NULL_CHAR) ? .25 : (*this)(nuc_hash(first_char), nuc_hash(*tmp1));
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


        /**
        * \brief Probability of the given amino acid sequence.
        */
//        prob_t aaProbability(const string& sequence, seq_len_t left_shift = 0, seq_len_t right_shift = 0) const {
//            return 0;
//        }


//        const event_matrix_t& nuc_prob_table() const {
//
//        }


//        const event_matrix_t& aa_prob_table() const {
//
//        }


//        const event_matrix_t& generateAminoAcidProbTable() {
//
//        }


//        bool readAminoAcidProbTable() {
//            return false;
//        }
//
//
//        bool writeAminoAcidProbTable() const {
//            return false;
//        }


        /**
        * \name Update chain probabilities.
        *
        * \brief Update internal vector of nucleotide joint probabilities.
        *
        * \param start Iterator to vector, in which each 4 elements are rows in nucleotide joint probability matrix, i.e.
        * row for prev-A-nucleotide, row for prev-C-nucleotide and so on.
        * \param mat Matrix with nucleotide joint probabilities with rows for prev nucleotide and columns
        * for next nucleotide.
        */
        ///@{
        void updateProbabilities(vector<prob_t>::const_iterator start) {
            if (_type == MonoNucleotide) {
                for (uint8_t i = 0; i < 4; ++i, ++start) {
                    this->_arr[i] = *start;
                }
            } else {
                for (uint8_t i = 0; i < 16; ++i, ++start) {
                    this->_arr[i] = *start;
                }
            }
        }

        void updateProbabilities(const event_matrix_t& mat) {
#ifdef YDEBUG
            if (_type == MonoNucleotide) { throw(std::runtime_error("Can't initialise the mono-nucleotide insertion model with a matrix!")); }
#endif

            for (uint8_t i = 0; i < 4; ++i) {
                for (uint8_t j = 0; j < 4; ++j) {
                    this->_arr[4*i + j] = mat(i, j);
                }
            }
        }
        ///@}


        /**
        * \name Event probability access.
        */
        ///@{
        prob_t operator[](uint8_t index) const { return this->_arr[index]; }
        prob_t operator()(uint8_t row, uint8_t col) const { return this->_arr[4*row + col]; }
        ///@}


        // update parameters given the AssemblyGraphRepertoire / AsemblyScenarioMatrix / iterator
        // getNucleotideSequenceProbability
        // getAminoAcidSequenceProbability
        // generate amino acid table
        // read amino acid table from file
        // read markov model w/ nucleotides from prob table
        // initialise markov model with input prob table / matrix
        // write ... to prob table (?)

    protected:

        INSERTION_MODEL_TYPE _type;
        prob_t* _arr;  /** Representation for a vector (#elements = 4) or a matrix (#elements = 16) with transition probabilities; rows and cols are for A-C-G-T (sequentially). */


        InsertionModel() {
            this->_arr = new prob_t[4];
            _type = MonoNucleotide;
        }

    };
}

#endif //YMIR_INSERTIONMODEL_H
