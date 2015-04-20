/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vdn at mailbox dot com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _MARKOV_CHAIN_H_
#define _MARKOV_CHAIN_H_


#include "string"

#include "types.h"


namespace ymir {

    /**
    * \class MarkovChain
    *
    * \brief Class for Markov Chain representing generation of sequences with joint distribution of neighbor symbols.
    */
    class MarkovChain {
    public:

        MarkovChain() {
            this->_arr = new prob_t[16];
        }


        MarkovChain(vector<prob_t>::const_iterator start) : MarkovChain() {
            this->updateProbabilities(start);
        }


        MarkovChain(const event_matrix_t& mat) : MarkovChain() {
            this->updateProbabilities(mat);
        }


        ~MarkovChain() {
            delete [] this->_arr;
        }


        /**
        * \brief Probability of the given nucleotide sequence.
        */
//        prob_t nucProbability(const string& sequence, char zero_symbol) const {
//            prob_t res = (*this)(nuc_hash(zero_symbol), nuc_hash(sequence[0]));
//            for (seq_len_t i = 1; i < sequence.size(); ++i) {
//                res *= (*this)(nuc_hash(sequence[i - 1]), nuc_hash(sequence[i]));
//            }
//            return res;
//        }


        /**
        * \brief Probability of the given nucleotide sequence.
        */
        prob_t nucProbability(const string& sequence) const {
            prob_t res = 1;
            for (seq_len_t i = 1; i < sequence.size(); ++i) {
                res *= (*this)(nuc_hash(sequence[i - 1]), nuc_hash(sequence[i]));
            }
//            return sequence.size() ? res : 0;
            return res;
        }


        prob_t nucProbability(string::const_iterator start, seq_len_t sequence_len) const {
            prob_t res = 1;
            string::const_iterator next = start + 1;
            for (seq_len_t i = 1; i < sequence_len; ++i, ++start, ++next) {
                res *= (*this)(nuc_hash(*start), nuc_hash(*next));
            }
//            return sequence_len ? res : 0;
            return res;
        }


        /**
        * \brief Probability of the given amino acid sequence.
        */
        prob_t aaProbability(const string& sequence, seq_len_t left_shift = 0, seq_len_t right_shift = 0) const {
            return 0;
        }


//        const event_matrix_t& nuc_prob_table() const {
//
//        }


//        const event_matrix_t& aa_prob_table() const {
//
//        }


//        const event_matrix_t& generateAminoAcidProbTable() {
//
//        }


        bool readAminoAcidProbTable() {
            return false;
        }


        bool writeAminoAcidProbTable() const {
            return false;
        }


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
            for (uint8_t i = 0; i < 16; ++i, ++start) {
                this->_arr[i] = *start;
            }
        }
        void updateProbabilities(const event_matrix_t& mat) {
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

        prob_t* _arr;  /** Representation for matrix (#elements = 16) with transition probabilities; rows and cols are for A-C-G-T (sequentially). */

    };
}

#endif