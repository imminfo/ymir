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


#include <unordered_map>
#include <iostream>
#include "tuple"
#include "matrix.h"

//#include "Eigen/Dense"

//#include "tools.h"

//#define MPFR
#include "gmp.h"
#include "mpfr.h"
//#include <Eigen/unsupported/Eigen/MPRealSupport>
//using namespace mpfr;
//using namespace Eigen;


namespace ymir {

    #define DEFAULT_DIV_GENE_MIN_LEN 2

    #define YDEBUG

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
    typedef uint16_t eventind_t;


    typedef uint16_t seq_len_t;


    /**
    * \typedef event_matrix_t
    *
    * \brief Matrix of stored event probabilities.
    */
//    typedef Eigen::Matrix<prob_t, Eigen::Dynamic, Eigen::Dynamic> event_matrix_t;
    typedef Matrix<prob_t, seq_len_t> event_matrix_t;


    typedef uint8_t segindex_t;


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


    typedef std::pair<eventind_t, prob_t> event_pair_t;


    /**
     * \enum MAAG_COMPUTE_PROB_ACTION
     */
    enum MAAG_COMPUTE_PROB_ACTION {
        MAX_PROBABILITY,
        SUM_PROBABILITY
    };


    enum GENE_SEGMENTS {
        VARIABLE,
        JOINING,
        DIVERSITY
    };


    /**
     * \enum MODEL_VECTOR_TYPE
     */
    enum RECOMBINATION {
        VJ_RECOMB,
        VDJ_RECOMB
    };


    /**
     * \enum EVENT_CLASS
     */
    enum EVENT_CLASS {
        NULL_EVENT = 0,

        VJ_VAR_JOI_GEN = 1,
        VJ_VAR_DEL = 2,
        VJ_JOI_DEL = 3,
        VJ_VAR_JOI_INS_LEN = 4,
        VJ_VAR_JOI_INS_NUC = 5,
        VJ_VAR_JOI_INS_NUC_A_ROW = 5,
        VJ_VAR_JOI_INS_NUC_C_ROW = 6,
        VJ_VAR_JOI_INS_NUC_G_ROW = 7,
        VJ_VAR_JOI_INS_NUC_T_ROW = 8,
        VJ_HYPMUT = 9,  // <- REMOVE THAT!!

        VDJ_VAR_GEN = 1,
        VDJ_JOI_DIV_GEN = 2,
        VDJ_VAR_DEL = 3,
        VDJ_JOI_DEL = 4,
        VDJ_DIV_DEL = 5,
        VDJ_VAR_DIV_INS_LEN = 6,
        VDJ_DIV_JOI_INS_LEN = 7,
        VDJ_VAR_DIV_INS_NUC = 8,
        VDJ_VAR_DIV_INS_NUC_A_ROW = 8,
        VDJ_VAR_DIV_INS_NUC_C_ROW = 9,
        VDJ_VAR_DIV_INS_NUC_G_ROW = 10,
        VDJ_VAR_DIV_INS_NUC_T_ROW = 11,
        VDJ_DIV_JOI_INS_NUC = 12,
        VDJ_DIV_JOI_INS_NUC_A_ROW = 12,
        VDJ_DIV_JOI_INS_NUC_C_ROW = 13,
        VDJ_DIV_JOI_INS_NUC_G_ROW = 14,
        VDJ_DIV_JOI_INS_NUC_T_ROW = 15,
        VDJ_HYPMUT = 16 // <- REMOVE THAT!!
    };


    enum MAAG_NODE_EVENT_INDEX {
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
}

#endif