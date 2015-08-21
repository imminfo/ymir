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

//    #ifndef DNDEBUG
    #define YDEBUG
//    #endif

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
     * \brief 1-based alignemnt of a Diversity gene to a sequence.
     */
    struct d_alignment_t {
        seq_len_t Dstart, Dend, seqstart, seqend;

        d_alignment_t() :
                Dstart(0), Dend(0), seqstart(0), seqend(0) {}


        d_alignment_t(seq_len_t Dstart, seq_len_t Dend, seq_len_t seqstart, seq_len_t seqend) :
                Dstart(Dstart), Dend(Dend), seqstart(seqstart), seqend(seqend) {}


        d_alignment_t(seq_len_t *p) :
                Dstart(*p), Dend(*(p + 1)), seqstart(*(p + 2)), seqend(*(p + 3)) {}


        bool operator==(const d_alignment_t& other) {
            return this->Dstart == other.Dstart
                    && this->Dend == other.Dend
                    && this->seqstart == other.seqstart
                    && this->seqend == other.seqend;
        }


        bool operator!=(const d_alignment_t& other) {
            return !((*this) == other);
        }

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
        UNKNOWN,
        VARIABLE,
        JOINING,
        DIVERSITY
    };


    /**
     * \enum Recombination
     */
    enum Recombination {
        UNDEFINED,
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


//    struct AlignmentsWithErrors {
//    public:
//
//        /**
//         * \struct vj_err_alignment_t
//         */
//        struct vj_err_alignment_t {
//            seq_len_t index;
//            event_ind_t event_index;
//            seq_len_t n_errors;
//
//
//            vj_err_alignment_t(seq_len_t index, event_ind_t event_index, seq_len_t n_errors)
//                    : index(index), event_index(event_index), n_errors(n_errors)
//            { }
//
//        };
//
//
//        /**
//         * \struct d_err_alignment_t
//         */
//        struct d_err_alignment_t {
//            seq_len_t row;
//            seq_len_t column;
//            event_ind_t event_index;
//            seq_len_t n_errors;
//
//
//            d_err_alignment_t(seq_len_t row, seq_len_t column, event_ind_t event_index, seq_len_t n_errors)
//                    : row(row), column(column), event_index(event_index), n_errors(n_errors)
//            { }
//
//        };
//
//
//        /**
//         * \brief Default constructor, nothing interesting.
//         */
//        AlignmentsWithErrors() { }
//
//
//        /**
//         * \brief Constructor which reserves memory for storing future alignments.
//         */
//        AlignmentsWithErrors(Recombination recomb) {
//            _vvec.reserve(DEFAULT_AWE_V_RESERVE_SIZE);
//            _jvec.reserve(DEFAULT_AWE_J_RESERVE_SIZE);
//            if (recomb == VDJ_RECOMB) {
//                _dvec.reserve(DEFAULT_AWE_D_RESERVE_SIZE);
//            }
//        }
//
//
//        /**
//         * \brief Add new alignment with errors to the storage.
//         *
//         */
//        ///@{
//        void addV(seq_len_t index, event_ind_t event_index, seq_len_t errors_num) {
//            _vvec.push_back(vj_err_alignment_t(index, event_index, errors_num));
//        }
//
//        void addJ(seq_len_t index, event_ind_t event_index, seq_len_t errors_num) {
//            _jvec.push_back(vj_err_alignment_t(index, event_index, errors_num));
//        }
//
//        void addD(seq_len_t row, seq_len_t column, event_ind_t event_index, seq_len_t errors_num) {
//            _dvec.push_back(d_err_alignment_t(row, column, event_index, errors_num));
//        }
//        ///@}
//
//
//        void finish() {
//            _vvec.resize(_vvec.size());
//            _jvec.resize(_jvec.size());
//            _dvec.resize(_dvec.size());
//        }
//
//
//        /**
//         * \brief Access to vectors with alignments with errors.
//         *
//         * \param i Index of the alignment.
//         */
//        ///@{
//        const vj_err_alignment_t& getV(size_t i) const { return _vvec[i]; }
//        const vj_err_alignment_t& getJ(size_t i) const { return _jvec[i]; }
//        const d_err_alignment_t& getD(size_t i) const { return _dvec[i]; }
//        ///@}
//
//
//        /**
//         * \brief Size of vectors with alignments.
//         */
//        ///@{
//        size_t sizeV() const { return _vvec.size(); }
//        size_t sizeJ() const { return _jvec.size(); }
//        size_t sizeD() const { return _dvec.size(); }
//        ///@}
//
//
//    protected:
//
//        std::vector<vj_err_alignment_t> _vvec;
//        std::vector<vj_err_alignment_t> _jvec;
//        std::vector<d_err_alignment_t> _dvec;
//
//    };


    typedef std::pair<std::string*, uint> codons_t;
    codons_t codons(char aminoacid) {
        switch (aminoacid) {
            default: return codons_t(nullptr, 0);
        }
    }
}

#endif