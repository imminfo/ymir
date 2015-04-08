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

#ifndef _MODELPARAMETERVECTOR_H_
#define _MODELPARAMETERVECTOR_H_

#define DEFAULT_DIV_GENE_MIN_LEN 3


#include <vector>

#include "types.h"


using namespace std;


namespace ymir {

    /**
     * \enum MODEL_VECTOR_TYPE
     */
    enum MODEL_VECTOR_TYPE {
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
        VJ_HYPMUT = 9,

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
        VDJ_HYPMUT = 16
    };


    /**
    * \class ModelParameterVector
    *
    * \brief Class for storing parameters of assembling statistical model. Note:
    * event with index 0 (zero) is "null" event and always has zero probability. Vector is ordered
    * so for VJ or VDJ recombination specific family events stored in specific order.
    *
    * Vector with event's probabilities stored in specific order:
    * VJ - recombination:
    *   [0] - null family
    *   [1] - Variable gene segment probability
    *   [2] - J prob
    *   [3 : 3 + [2] - [1]] - V deletions
    *   [4] - J deletions
    *   [5] - VJ insertion length
    *   [6] - VJ insertion markov chain
    *
    * VDJ - recombination:
    *   [0] - null family
    *   [1] - Variable gene segment probability
    *   [2] - J prob
    *   [3] - V deletions
    *   [4] - J deletions
    *   [5] - VJ insertion length
    *   [6] - VJ insertion markov chain
    *
    * Hypermutations added as the last element in the vector.
    * Note: gene segment deletion probabilities stored in order of gene segment probabilities.
    */
    class ModelParameterVector {
    public:


        /**
         * Event family - family of specific events like deletions of V4 segment or DJ insertions.
         * Event class - class of events like V deletions, J deletions, D deletions, VJ insertions.
         *
         * \param param_vec Vector of probabilities of events.
         * \param lens_vec Vector of lengths for each family of events (V1 deletions, V2 deletions, etc.)
         * \param event_classes Vector of start indices in lens_vec of each class of events (all V deletions, all J deletions, VJ insertions, etc.)
         * \param event_family_row_numbers Number of rows in each event family (zero means there is no matrix, but only a vector). Vector with the same
         * length as lens_vec.
         * \param laplace_vec Vector of pseudo counts for each event family.
         */
        ModelParameterVector(MODEL_VECTOR_TYPE vec_type,
                             const vector<prob_t>& param_vec,
                             const vector<eventind_t>& lens_vec,
                             const vector<eventind_t>& event_classes,
                             const vector<seq_len_t>& event_family_row_numbers,
                             const vector<prob_t>& laplace_vec = vector<prob_t>(),
                             bool do_normalise = true,
                             const vector<seq_len_t>& d_genes_min_len = vector<seq_len_t>()) {
            _vec = vector<prob_t>();
            _vec.reserve(param_vec.size() + 1);
            _vec.push_back(0);
            _vec.insert(_vec.end(), param_vec.begin(), param_vec.end());

            _event_classes = vector<eventind_t>();
            _event_classes.reserve(event_classes.size() + 1);
            for (eventind_t i = 0; i < event_classes.size(); ++i) {
                _event_classes.push_back(event_classes[i] + 1);
            }

            _event_family_row_numbers = vector<seq_len_t>();
            _event_family_row_numbers.reserve(event_family_row_numbers.size());
            _event_family_row_numbers.insert(_event_family_row_numbers.end(),
                                            event_family_row_numbers.begin(),
                                            event_family_row_numbers.end());

            _edges.reserve(lens_vec.size() + 4);
            _edges.push_back(0);
            _edges.push_back(1);
            for (eventind_t i = 0; i < lens_vec.size(); ++i) {
                _edges.push_back(_edges[i + 1] + lens_vec[i]);
            }
            _edges.push_back(_vec.size());
            _event_classes.push_back(_edges.size() - 1);

            _laplace = vector<prob_t>();
            if (laplace_vec.size()) {
                _laplace.reserve(laplace_vec.size());
                _laplace.push_back(0);
                for (eventind_t i = 0; i < laplace_vec.size(); ++i) {
                    _laplace.push_back(laplace_vec[i]);
                }
            } else {
                // Laplace correction equal to zero if vector is not supplied.
                _laplace.resize(_edges.size(), 0);
            }

            if (vec_type == VDJ_RECOMB) {
                _d_genes_min_len = vector<seq_len_t>();
                if (d_genes_min_len.size()) {
                    _d_genes_min_len = d_genes_min_len;
                } else {
                    _d_genes_min_len.resize(this->eventClassSize(VDJ_DIV_DEL), DEFAULT_DIV_GENE_MIN_LEN);
                }
            }

            if (do_normalise) {
                this->normaliseEventFamilies();
            }
        }


        //============= VECTOR INDICES ACCESS =============//


        eventind_t eventFamilySize(eventind_t event_family) const {
            return _edges[event_family + 1] - _edges[event_family];
        }

        eventind_t eventFamilySize(eventind_t event_class, eventind_t event_family) const {
            return _edges[_event_classes[event_class] + 1] - _edges[_event_classes[event_class]];
        }

        eventind_t eventClassSize(eventind_t event_class) const {
            return _edges[_event_classes[event_class + 1]] - _edges[_event_classes[event_class]];
        }


        /**
         * \brief Get a probability of an event with the given global index.
         *
         * \param event_family Event family's index.
         * \param gl_event_index Global index of the event.
         * \param loc_event_index Local index of the event, i.e., index of the event in the given event family.
         *
         * \return Probability of the event.
         */
        prob_t operator[] (eventind_t gl_event_index) const { return _vec[gl_event_index]; }

        prob_t operator() (eventind_t event_family, eventind_t loc_event_index) const {
            return _vec[_edges[event_family] + loc_event_index];
        }

        vector<prob_t>::const_iterator get_iterator(eventind_t i) const { return _vec.begin() + i; }


        /**
         *
         * // vector of events
         * param_vec(VDJ_VAR_GEN, 5) // 5th V segment
         * vec_ind(VDJ_VAR_GEN, 5)
         *
         * // two different event classes with similar access
         * // matrix of events
         * param_vec(VJ_VAR_JOI_GEN, 3, 4)  // 3th V segment, 4th J segment
         * mat_ind(VJ_VAR_JOI_GEN, 3, 4)
         * // vector of vectors of events
         * param_vec(VJ_VAR_DEL, 5, 15)  // 5th V segment, 15 deletions
         * vec_ind(VJ_VAR_DEL, 5, 15)
         * // also
         * // matrix of events
         * param_vec(VDJ_DIV_DEL, 2, 7, 8) // 2nd D segment, 7 D3' deletions, 8 D5' deletions
         * mat_ind(VDJ_DIV_DEL, 2, 7, 8)
         */
        // get indices from vector of events
        eventind_t event_index(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_index) const {
            return _edges[_event_classes[event_class] + (event_family - 1)] + event_index;
        }
        eventind_t event_index(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_row, eventind_t event_column) const {
            return _edges[_event_classes[event_class] + (event_family - 1)]
                   + (event_row - 1) * _event_family_row_numbers[_event_classes[event_class] + (event_family - 1)] + (event_column - 1);
        }
        eventind_t vec_index(EVENT_CLASS event_class, eventind_t event_index) const {
            return _edges[_event_classes[event_class]] + event_index;
        }
        // get indices from specific vector of events from ordered set of vectors
        eventind_t vec_index(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_index) const {
            return _edges[_event_classes[event_class] + (event_family - 1)] + event_index;
        }
        // get indices from matrix of events
        eventind_t mat_index(EVENT_CLASS event_class, eventind_t event_row, eventind_t event_column) const {
            return _edges[_event_classes[event_class]]
                   + (event_row - 1) * _event_family_row_numbers[_event_classes[event_class]] + (event_column - 1);
        }
        // get indices specific matrix from ordereded set of matrices of events
        eventind_t mat_index(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_row, eventind_t event_column) const {
            return _edges[_event_classes[event_class] + (event_family - 1)]
                   + (event_row - 1) * _event_family_row_numbers[_event_classes[event_class] + (event_family - 1)] + (event_column - 1);
        }

        prob_t event_prob(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_index) const {
            return _vec[this->event_index(event_class, event_family, event_index)];
        }
        prob_t event_prob(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_row, eventind_t event_column) const {
            return _vec[this->event_index(event_class, event_family, event_row, event_column)];
        }
        prob_t vec_prob(EVENT_CLASS event_class, eventind_t event_index) const {
            return (*this)[vec_index(event_class, event_index)];
        }
        prob_t vec_prob(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_index) const {
            return (*this)[vec_index(event_class, event_family, event_index)];
        }
        prob_t mat_prob(EVENT_CLASS event_class, eventind_t event_row, eventind_t event_column) const {
            return (*this)[mat_index(event_class, event_row, event_column)];
        }
        prob_t mat_prob(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_row, eventind_t event_column) const {
            return (*this)[mat_index(event_class, event_family, event_row, event_column)];
        }


        /**
        * \brief Normalise each event family to have sum equal to 1.
        */
        void normaliseEventFamilies() {
            for (eventind_t i = 2; i < _edges.size(); ++i) {
                prob_t prob_sum = _laplace[i-1] * (_edges[i] - _edges[i-1]);
                for (eventind_t j = _edges[i-1]; j < _edges[i]; ++j) {
                    prob_sum += _vec[j];
                }
                for (eventind_t j = _edges[i-1]; j < _edges[i]; ++j) {
                    _vec[j] = (_vec[j] + _laplace[i-1]) / prob_sum;
                }
            }
        }


        /**
         *
         */
        eventind_t max_VJ_ins_len() const {
            return this->eventFamilySize(VJ_VAR_JOI_INS_LEN, 0) - 1;
        }
        eventind_t max_VD_ins_len() const {
            return this->eventFamilySize(VDJ_VAR_DIV_INS_LEN, 0) - 1;
        }
        eventind_t max_DJ_ins_len() const {
            return this->eventFamilySize(VDJ_DIV_JOI_INS_LEN, 0) - 1;
        }


        /**
         *
         */
        inline seq_len_t D_min_len(eventind_t d_index) const { return _d_genes_min_len[d_index - 1]; }

    private:

        vector<prob_t> _vec;
        vector<eventind_t> _edges;  /** Vector with starting indices for each event family. */
        vector<eventind_t> _event_classes;  /** Vector of indices of event classes in _edges. */
        vector<seq_len_t> _event_family_row_numbers;  /** Vector of the number of rows of each event family. */
        vector<prob_t> _laplace;
        vector<seq_len_t> _d_genes_min_len;


        /**
        * \brief Private default constructor.
        */
        ModelParameterVector() {}

    };
}

#endif