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
         * \param event_family_col_numbers Number of rows in each event family (zero means there is no matrix, but only a vector). Vector with the same
         * length as lens_vec.
         * \param laplace_vec Vector of pseudo counts for each event family.
         */
        ModelParameterVector(MODEL_VECTOR_TYPE vec_type,
                             const vector<prob_t>& param_vec,
                             const vector<eventind_t>& lens_vec,
                             const vector<eventind_t>& event_classes,
                             const vector<seq_len_t>& event_family_col_numbers,
                             const vector<prob_t>& laplace_vec = vector<prob_t>(),
                             bool do_normalise = true,
                             const vector<seq_len_t>& d_genes_min_len = vector<seq_len_t>()) {
            _vec = vector<prob_t>();
            _vec.reserve(param_vec.size() + 1);
            _vec.push_back(0);
            _vec.insert(_vec.end(), param_vec.begin(), param_vec.end());

            _event_classes = vector<eventind_t>();
            _event_classes.reserve(event_classes.size() + 1);
            _event_classes.push_back(0);
            for (eventind_t i = 0; i < event_classes.size(); ++i) {
                _event_classes.push_back(event_classes[i] + 1);
            }

            _event_family_col_numbers = vector<seq_len_t>();
            _event_family_col_numbers.reserve(event_family_col_numbers.size() + 1);
            _event_family_col_numbers.push_back(0);
            for (eventind_t i = 0; i < event_family_col_numbers.size(); ++i) {
                _event_family_col_numbers.push_back(event_family_col_numbers[i]);
            }

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


        bool operator==(const ModelParameterVector &other) const {
            return (_vec == other._vec)
                   && (_edges == other._edges)
                   && (_event_classes == other._event_classes)
                   && (_event_family_col_numbers == other._event_family_col_numbers)
                   && (_laplace == other._laplace)
                   && (_d_genes_min_len == other._d_genes_min_len);
        }


        //============= VECTOR INDICES ACCESS =============//


        eventind_t eventFamilySize(eventind_t event_family) const {
            return _edges[event_family + 1] - _edges[event_family];
        }

        eventind_t eventFamilySize(eventind_t event_class, eventind_t event_family) const {
            return eventFamilySize(_event_classes[event_class] + event_family);
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
         * 0-based indices for families and events! All segment indices must be passed as (segindex - 1),
         * segment families as (segindex_deletion_index - 1), but deletions and insertions should be passed as it is (#deletions).
         */
        eventind_t event_index(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_index) const {
            return _edges[_event_classes[event_class] + event_family] + event_index;
        }
        eventind_t event_index(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_row, eventind_t event_column) const {
            return _edges[_event_classes[event_class] + event_family]
                   + event_row * _event_family_col_numbers[_event_classes[event_class] + event_family] + event_column;
        }

        prob_t event_prob(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_index) const {
            return _vec[this->event_index(event_class, event_family, event_index)];
        }
        prob_t event_prob(EVENT_CLASS event_class, eventind_t event_family, eventind_t event_row, eventind_t event_column) const {
            return _vec[this->event_index(event_class, event_family, event_row, event_column)];
        }



        //============= EVENT ACCESS =============//



//        eventind_t index_V_gene(segindex_t v_index) const { return event_index(); }
//        prob_t prob_V_gene(segindex_t v_index) const { return _vec[index_V_gene(v_index)]; }  // Hmmm...
//
//
//        eventind_t index_VJ_genes(segindex_t v_index, segindex_t j_index) const {
//            return _edges[1] + (v_index - 1) * _v_gene_num + (j_index - 1);
//        }
//        prob_t prob_VJ_genes(segindex_t v_index, segindex_t j_index) const {
//            return _vec[index_VJ_genes(v_index, j_index)];
//        }  // Hmmm...
//
//
//        eventind_t index_JD_genes(segindex_t j_index, segindex_t d_index) const { return _edges[2] + (j_index - 1) * _d_gene_num + (d_index - 1); }
//        prob_t prob_JD_genes(segindex_t j_index, segindex_t d_index) const { return _vec[index_JD_genes(j_index, d_index)]; }
//
//
//        eventind_t index_V_del(segindex_t v_index, seq_len_t del_num) const { return _edges[3 + (v_index - 1)] + del_num; }
//        prob_t prob_V_del(segindex_t v_index, seq_len_t del_num) const { return _vec[index_V_del(v_index, del_num)]; }
//
//
//        eventind_t index_J_del(segindex_t j_index, seq_len_t del_num) const { return _edges[3 + (_edges[2] - 1) + (j_index - 1)] + del_num; }
//        prob_t prob_J_del(segindex_t j_index, seq_len_t del_num) const { return _vec[index_J_del(j_index, del_num)]; }
//
//
//        eventind_t index_D_del(segindex_t d_index, seq_len_t d5_del_num, seq_len_t d3_del_num) const {
//            return _edges[3 + (_edges[2] - 1) + 2 + (d_index - 1)] + d5_del_num*_d_gene_max_dels[d_index - 1] + d3_del_num;
//        }
//        prob_t prob_D_del(segindex_t d_index, seq_len_t d5_del_num, seq_len_t d3_del_num) const {
//            return _vec[index_D_del(d_index, d5_del_num, d3_del_num)];
//        }
//
//
//        eventind_t index_VJ_ins_len(seq_len_t ins_len) const {
//            return _edges[_edges.size() - 7] + ins_len; // - 1 - 4 - 2
//        }
//        prob_t prob_VJ_ins_len(seq_len_t ins_len) const {
//            return _vec[index_VJ_ins_len(ins_len)];
//        }
//
//
//        eventind_t index_VD_ins_len(seq_len_t ins_len) const {
//            return _edges[_edges.size() - 12] + ins_len; // - 1 - 8 - 3
//        }
//        prob_t prob_VD_ins_len(seq_len_t ins_len) const {
//            return _vec[index_VD_ins_len(ins_len)];
//        }
//
//
//        eventind_t index_DJ_ins_len(seq_len_t ins_len) const {
//            return _edges[_edges.size() - 11] + ins_len; // -1 - 8 - 2
//        }
//        prob_t prob_DJ_ins_len(seq_len_t ins_len) const {
//            return _vec[index_DJ_ins_len(ins_len)];
//        }


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
        vector<seq_len_t> _event_family_col_numbers;  /** Vector of the number of rows of each event family. */
        vector<prob_t> _laplace;
        vector<seq_len_t> _d_genes_min_len;


        /**
        * \brief Private default constructor.
        */
        ModelParameterVector() {}

    };
}

#endif