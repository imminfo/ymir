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
#include "tuple"

#include "Eigen/Dense"

#include "tools.h"


namespace ymir {
    /**
    * \typedef numeric
    *
    * \brief Type of assembly scenario probabilities.
    */
//    #define MPFR
    #ifdef MPFR
    #include "Eigen/mpreal.h"
    typedef mpfr::mpreal numeric;

    /**
    * \typedef prob_t
    *
    * \brief Type of stored probabilities of different events.
    */
    typedef mpfr::mpreal prob_t;
    #else
    /**
    * \typedef numeric
    *
    * \brief Type of sequences assembling probability.
    */
    typedef double numeric;

    /**
    * \typedef prob_t
    *
    * \brief Type of stored probabilities of different events.
    */
    typedef double prob_t;
    #endif


    /**
    * \typedef eventind_t
    *
    * \brief Type of stored indices of events in vertices. Zero means no event or some complex event.
    */
    typedef uint16_t eventind_t;


    /**
    * \typedef event_matrix_t
    *
    * \brief Matrix of stored event probabilities.
    */
    typedef Eigen::Matrix<prob_t, Eigen::Dynamic, Eigen::Dynamic> event_matrix_t;


    typedef uint16_t seq_len_t;


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


    struct NamedVectorArray {

    public:

        struct Vector {

        public:

            vector<prob_t> vec;
            string name;

            Vector(const string& name_) {
                name = name_;
                vec.clear();
                vec.reserve(10);
            }
        };

        NamedVectorArray() {
            cols.clear();
            cols.reserve(10);
        }

        void addColumn(const string& name) {
            cols.push_back(name);
            map[name] = cols.size() - 1;
        }

        Vector& operator[](const string& name) {
            return cols[map[name]];
        }

        void push(int index, prob_t value) {
            cols[index].vec.push_back(value);
        }

        string getName(int index) const {
            return cols[index].name;
        }


        vector<Vector> cols;
        unordered_map<string, int> map;
    };
}

#endif