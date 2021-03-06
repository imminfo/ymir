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
#include <vector>

//#include "jsoncpp.cpp"
#include "json.hpp"
using json = nlohmann::json;

#ifdef USE_OMP
#include "omp.h"
#endif

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
//    #define YDEBUG
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
//    typedef double prob_t;
    typedef long double prob_t;
    #endif


    typedef std::string sequence_t;


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


    using std::unique_ptr;


    using std::shared_ptr;


    typedef std::pair<event_ind_t, prob_t> event_pair_t;


    typedef uint8_t codon_hash;

//    typedef Json::Value json_value;  // jsoncpp.h
    typedef json json_value;  // json.hpp

    json_value get_json_value(const json_value &obj, const json::object_t::key_type &key, json::value_t default_value) {
        // jsoncpp.cpp
        // return obj.get(key, default_value);
        // json.hpp
        return obj.value(key, default_value);
    }

    json_value get_json_value(const json_value &obj, const json::object_t::key_type &key, const char* default_value) {
        return json_value(obj.value(key, default_value));
    }

    json_value get_json_value(const json_value &obj, const json::object_t::key_type &key, bool default_value) {
        return obj.value(key, default_value);
    }

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
        UNDEF_SEQ_TYPE = 0,
        NUCLEOTIDE = 1,
        AMINOACID = 2
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


    enum SequenceCodingType {
        ALL = 0,
        CODING = 1,
        NONCODING = 2,
        OUTOFFRAME = 3,
        STOPCODON = 4
    };

}

#endif