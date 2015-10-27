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

#ifndef _TYPES_H
#define _TYPES_H


#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "jsoncpp.cpp"

#include "matrix.h"

//#include "Eigen/Dense"

//#include "tools.h"

//#define MPFR
//#include "gmp.h"
//#include "mpfr.h"
//#include <Eigen/unsupported/Eigen/MPRealSupport>
//using namespace mpfr;
//using namespace Eigen;


namespace ymir {

    #define DEFAULT_DIV_GENE_MIN_LEN 3

    #define NULL_CHAR '_'

    #define DEFAULT_MAX_INS_LENGTH 65

    #ifndef DNDEBUG
    #define YDEBUG
    #endif

    #define DEFAULT_AWE_V_RESERVE_SIZE 60
    #define DEFAULT_AWE_D_RESERVE_SIZE 4000
    #define DEFAULT_AWE_J_RESERVE_SIZE 60

    /**
    * \typedef numeric
    *
    * \brief Type of assembly scenario probabilities.
    */

    #ifdef MPFR

    /**
    * \typedef prob_t
    *
    * \brief Type of stored probabilities of different events.
    */
    typedef mpreal prob_t;
    #else
    /**
    * \typedef prob_t
    *
    * \brief Type of stored probabilities of different events.
    */
    typedef double prob_t;
//    typedef long double prob_t;
    #endif


    /**
    * \typedef eventind_t
    *
    * \brief Type of stored indices of events in vertices. Zero means no event or some complex event.
    */
    typedef uint16_t event_ind_t;


    typedef uint16_t seq_len_t;


    /**
    * \typedef event_matrix_t
    *
    * \brief Matrix of stored event probabilities.
    */
//    typedef Eigen::Matrix<prob_t, Eigen::Dynamic, Eigen::Dynamic> event_matrix_t;
    typedef Matrix<prob_t, seq_len_t> event_matrix_t;


    typedef uint8_t seg_index_t;


    /**
     * \struct d_alignment_t
     *
     * \brief 1-based alignment of a gene segment sequence to a reference sequence without indels.
     */
    class Alignment {

    public:

        Alignment()
                : _gene_start(0), _seq_start(0), _len(0)
        { }


        Alignment(seq_len_t gene_start, seq_len_t seq_start, seq_len_t alignment_length)
                : _gene_start(gene_start), _seq_start(seq_start), _len(alignment_length)
        { }


        Alignment(seq_len_t *p)
                : _gene_start(*p), _seq_start(*(p + 1)), _len(*(p + 2))
        { }


        ///@{
        seq_len_t gene_start() const { return _gene_start; }

        seq_len_t gene_end() const { return _gene_start + _len - 1; }

        seq_len_t seq_start() const { return _seq_start; }

        seq_len_t seq_end() const { return _seq_start + _len - 1; }

        seq_len_t length() const { return _len; }
        ///@}


        bool operator==(const Alignment& other) {
            return this->_gene_start == other._gene_start
                   && this->_seq_start == other._seq_start
                   && this->_len == other._len;
        }


        bool operator!=(const Alignment& other) {
            return !((*this) == other);
        }

    private:

        seq_len_t _gene_start, _seq_start, _len;

    };


    typedef std::pair<event_ind_t, prob_t> event_pair_t;


    /**
     * \enum MAAG_COMPUTE_PROB_ACTION
     */
    enum MAAGComputeProbAction {
        MAX_PROBABILITY,
        SUM_PROBABILITY
    };


    enum GeneSegments {
        UNDEF_GENE,
        VARIABLE,
        JOINING,
        DIVERSITY
    };


    /**
     * \enum Recombination
     */
    enum Recombination {
        UNDEF_RECOMB,
        VJ_RECOMB,
        VDJ_RECOMB
    };


    /**
     * \enum ModelBehaviour
     */
    enum ModelBehaviour {
        PREDEFINED,
        EMPTY
    };


    /**
     * \enum MetadataMode
     */
    enum MetadataMode {
        NO_METADATA = 0,
        SAVE_METADATA = 1
    };


    /**
     * \enum ErrorMode
     */
    enum ErrorMode {
        NO_ERRORS = 0,
        COMPUTE_ERRORS = 1
    };


    /**
     * \enum SequenceType
     */
    enum SequenceType {
        NUCLEOTIDE,
        AMINOACID
    };


    /**
     * \enum EVENT_CLASS
     */
    enum EventClass {
        NULL_EVENT = 0,

        VJ_VAR_JOI_GEN = 1,
        VJ_VAR_DEL = 2,
        VJ_JOI_DEL = 3,
        VJ_VAR_JOI_INS_LEN = 4,
        VJ_VAR_JOI_INS_NUC = 5,
        VJ_VAR_JOI_INS_NUC_A = 5,
        VJ_VAR_JOI_INS_NUC_C = 6,
        VJ_VAR_JOI_INS_NUC_G = 7,
        VJ_VAR_JOI_INS_NUC_T = 8,
        VJ_ERROR_RATE = 9,

        VDJ_VAR_GEN = 1,
        VDJ_JOI_DIV_GEN = 2,
        VDJ_VAR_DEL = 3,
        VDJ_JOI_DEL = 4,
        VDJ_DIV_DEL = 5,
        VDJ_VAR_DIV_INS_LEN = 6,
        VDJ_DIV_JOI_INS_LEN = 7,
        VDJ_VAR_DIV_INS_NUC = 8,
        VDJ_VAR_DIV_INS_NUC_A = 8,
        VDJ_VAR_DIV_INS_NUC_C = 9,
        VDJ_VAR_DIV_INS_NUC_G = 10,
        VDJ_VAR_DIV_INS_NUC_T = 11,
        VDJ_DIV_JOI_INS_NUC = 12,
        VDJ_DIV_JOI_INS_NUC_A = 12,
        VDJ_DIV_JOI_INS_NUC_C = 13,
        VDJ_DIV_JOI_INS_NUC_G = 14,
        VDJ_DIV_JOI_INS_NUC_T = 15,
        VDJ_ERROR_RATE = 16
    };


    enum MAAGNodeEventIndex {
        VJ_VAR_JOI_GEN_I = 0,
        VJ_VAR_DEL_I = 1,
        VJ_VAR_JOI_INS_I = 2,
        VJ_JOI_DEL_I = 3,

        VDJ_VAR_GEN_I = 0,
        VDJ_VAR_DEL_I = 1,
        VDJ_VAR_DIV_INS_I = 2,
        VDJ_DIV_DEL_I = 3,
        VDJ_DIV_JOI_INS_I = 4,
        VDJ_JOI_DEL_I = 5,
        VDJ_JOI_DIV_GEN_I = 6
    };


    enum InsertionModelType {
        MONO_NUCLEOTIDE,
        DI_NUCLEOTIDE
    };


    /**
     * \struct CodonTable
     *
     * \brief A struct for representing nucleotide codons for amino acids.
     */
    struct CodonTable {

        struct Codons {

            Codons(std::pair<std::unordered_multimap<char, std::string>::const_iterator, std::unordered_multimap<char, std::string>::const_iterator> it)
                    : _begin(it.first), _end(it.second), _current(it.first)
            {}


            std::string next() {
                std::string res = _current->second;
                ++_current;
                return res;
            }


            bool end() const { return _current == _end; }

        protected:
            std::unordered_multimap<char, std::string>::const_iterator _begin, _end, _current;

            Codons() {}

        };

        CodonTable() {
            _codons = {
                    {'A', "GCT"}, {'A', "GCC"}, {'A', "GCA"}, {'A', "GCG"},
                    {'L', "TTA"}, {'L', "TTG"}, {'L', "CTT"}, {'L', "CTC"}, {'L', "CTA"}, {'L', "CTG"},
                    {'R', "CGT"}, {'R', "CGC"}, {'R', "CGA"}, {'R', "CGG"}, {'R', "AGA"}, {'R', "AGG"},
                    {'K', "AAA"}, {'K', "AAG"},
                    {'N', "AAT"}, {'N', "AAC"},
                    {'M', "ATG"},
                    {'D', "GAT"}, {'D', "GAC"},
                    {'F', "TTT"}, {'F', "TTC"},
                    {'C', "TGT"}, {'C', "TGC"},
                    {'P', "CCT"}, {'P', "CCC"}, {'P', "CCA"}, {'P', "CCG"},
                    {'Q', "CAA"}, {'Q', "CAG"},
                    {'S', "TCT"}, {'S', "TCC"}, {'S', "TCA"}, {'S', "TCG"}, {'S', "AGT"}, {'S', "AGC"},
                    {'E', "GAA"}, {'E', "GAG"},
                    {'T', "ACT"}, {'T', "ACC"}, {'T', "ACA"}, {'T', "ACG"},
                    {'G', "GGT"}, {'G', "GGC"}, {'G', "GGA"}, {'G', "GGG"},
                    {'W', "TGG"},
                    {'H', "CAT"}, {'H', "CAC"},
                    {'Y', "TAT"}, {'Y', "TAC"},
                    {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"},
                    {'V', "GTT"}, {'V', "GTC"}, {'V', "GTA"}, {'V', "GTG"},
                    {'*', "TAA"}, {'*', "TGA"}, {'*', "TAG"}
            };
        }

        Codons codons(char aminoacid) const { return Codons(_codons.equal_range(aminoacid)); }



    protected:
        std::unordered_multimap<char, std::string> _codons;
    };


    typedef std::pair<std::string*, uint> codons_t;
    codons_t codons(char aminoacid) {
        switch (aminoacid) {
            default: return codons_t(nullptr, 0);
        }
    }


    struct AlignmentEventScore {
        float match, mismatch, ins, del;

        AlignmentEventScore()
                : match(1), mismatch(-1), ins(-1), del(-1)
        {}

        AlignmentEventScore(float match_, float mismatch_, float ins_, float del_)
                : match(match_), mismatch(mismatch_), ins(ins_), del(del_)
        {}
    };


    /**
     * \brief Events in alignments match, mismatch, indels represented as sets of bits.
     */
    ///@{
    typedef std::bitset<2> alignment_event_t;

    const alignment_event_t MATCH = alignment_event_t(0);

    const alignment_event_t MISMATCH = alignment_event_t(1);

    const alignment_event_t INSERTION = alignment_event_t(2);

    const alignment_event_t DELETION = alignment_event_t(3);
    ///@}


    /**
     * \class AlignmentEventVector
     */
    class AlignmentEventVector {
    public:

        typedef std::vector<bool> storage_t;


        AlignmentEventVector() {}


        void addEvent(const alignment_event_t &event) {
            _vec.push_back(event[0]);
            _vec.push_back(event[1]);
        }


        alignment_event_t operator[](size_t index) const {
            return alignment_event_t(2*_vec[index*2] + _vec[index*2 + 1]);
        }

    private:

        storage_t _vec;

    };

}

#endif