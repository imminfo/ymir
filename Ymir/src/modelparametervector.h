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
        ModelParameterVector(Recombination vec_type,
                             const vector<prob_t>& param_vec,
                             const vector<event_ind_t>& lens_vec,
                             const vector<event_ind_t>& event_classes,
                             const vector<seq_len_t>& event_family_col_numbers,
                             const vector<prob_t>& laplace_vec = vector<prob_t>(),
                             bool do_normalise = true,
                             const vector<seq_len_t>& d_genes_min_len = vector<seq_len_t>()) {
            _recomb = vec_type;

            _vec = vector<prob_t>();
            _vec.reserve(param_vec.size() + 1);
            _vec.push_back(0);
            _vec.insert(_vec.end(), param_vec.begin(), param_vec.end());

            _event_classes = vector<event_ind_t>();
            _event_classes.reserve(event_classes.size() + 1);
            _event_classes.push_back(0);
            for (event_ind_t i = 0; i < event_classes.size(); ++i) {
                _event_classes.push_back(event_classes[i] + 1);
            }

            _event_family_col_numbers = vector<seq_len_t>();
            _event_family_col_numbers.reserve(event_family_col_numbers.size() + 1);
            _event_family_col_numbers.push_back(0);
            for (event_ind_t i = 0; i < event_family_col_numbers.size(); ++i) {
                _event_family_col_numbers.push_back(event_family_col_numbers[i]);
            }

            _edges.reserve(lens_vec.size() + 4);
            _edges.push_back(0);
            _edges.push_back(1);
            for (event_ind_t i = 0; i < lens_vec.size(); ++i) {
                _edges.push_back(_edges[i + 1] + lens_vec[i]);
            }
            _edges.push_back(_vec.size());
            _event_classes.push_back(_edges.size() - 1);

            _laplace = vector<prob_t>();
            if (laplace_vec.size()) {
                _laplace.reserve(laplace_vec.size() + 1);
                _laplace.push_back(0);
                for (event_ind_t i = 0; i < laplace_vec.size(); ++i) {
                    _laplace.push_back(laplace_vec[i]);
                }
            } else {
                // Laplace correction equal to zero if vector is not supplied.
                _laplace.resize(_edges.size() - 2, 0);
            }

            _d_genes_min_len = vector<seq_len_t>();
            if (vec_type == VDJ_RECOMB) {
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
            if ((_recomb != other._recomb)
                || (_vec.size() != other._vec.size())
                || (_edges.size() != other._edges.size())
                || (_event_classes.size() != other._event_classes.size())
                || (_event_family_col_numbers.size() != other._event_family_col_numbers.size())
                || (_laplace.size() != other._laplace.size())
                || (_d_genes_min_len.size() != other._d_genes_min_len.size())) {
                return false;
            }

            for (size_t i = 0; i < _vec.size(); ++i) { if (_vec[i] != other._vec[i]) { return false; } }
            for (size_t i = 0; i < _edges.size(); ++i) { if (_edges[i] != other._edges[i]) { return false; } }
            for (size_t i = 0; i < _event_classes.size(); ++i) { if (_event_classes[i] != other._event_classes[i]) { return false; } }
            for (size_t i = 0; i < _event_family_col_numbers.size(); ++i) { if (_event_family_col_numbers[i] != other._event_family_col_numbers[i]) { return false; } }
            for (size_t i = 0; i < _laplace.size(); ++i) { if (_laplace[i] != other._laplace[i]) { return false; } }
            for (size_t i = 0; i < _d_genes_min_len.size(); ++i) { /* cout << (int) _d_genes_min_len[i] << ":" << (int) other._d_genes_min_len[i] << endl; */ if (_d_genes_min_len[i] != other._d_genes_min_len[i]) { return false; } }

            return true;
        }


        Recombination recombination() const { return _recomb; }


        //============= VECTOR INDICES ACCESS =============//


        event_ind_t eventFamilySize(event_ind_t event_family) const {
            return _edges[event_family + 1] - _edges[event_family];
        }

        event_ind_t eventFamilySize(event_ind_t event_class, event_ind_t event_family) const {
            return eventFamilySize(_event_classes[event_class] + event_family);
        }

        event_ind_t eventClassSize(event_ind_t event_class) const {
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
        ///@{
        prob_t& operator[] (event_ind_t gl_event_index) { return _vec[gl_event_index]; }

        prob_t operator[] (event_ind_t gl_event_index) const { return _vec[gl_event_index]; }

        prob_t operator() (event_ind_t event_family, event_ind_t loc_event_index) const {
            return _vec[_edges[event_family] + loc_event_index];
        }
        ///@}

        vector<prob_t>::const_iterator get_iterator(event_ind_t i) const { return _vec.begin() + i; }


        /**
         * 0-based indices for families and events! All segment indices must be passed as (segindex - 1),
         * segment families as (segindex_deletion_index - 1), but deletions and insertions should be passed as it is (#deletions).
         */
        event_ind_t event_index(EventClass event_class, event_ind_t event_family, event_ind_t event_index) const {
            if (_edges[_event_classes[event_class] + event_family] + event_index < _edges[_event_classes[event_class] + event_family + 1]) {
                return _edges[_event_classes[event_class] + event_family] + event_index;
            }
            return 0;

//            return _edges[_event_classes[event_class] + event_family] + event_index;
        }

        event_ind_t event_index(EventClass event_class, event_ind_t event_family, event_ind_t event_row, event_ind_t event_column) const {
            if (event_column < _event_family_col_numbers[_event_classes[event_class] + event_family] && _edges[_event_classes[event_class] + event_family]
                + event_row * _event_family_col_numbers[_event_classes[event_class] + event_family] + event_column < _edges[_event_classes[event_class] + event_family + 1]) {
                return _edges[_event_classes[event_class] + event_family]
                       + event_row * _event_family_col_numbers[_event_classes[event_class] + event_family] + event_column;
            }
            return 0;

//            return _edges[_event_classes[event_class] + event_family]
//                   + event_row * _event_family_col_numbers[_event_classes[event_class] + event_family] + event_column;
        }

        prob_t event_prob(EventClass event_class, event_ind_t event_family, event_ind_t event_index) const {
            return _vec[this->event_index(event_class, event_family, event_index)];
        }

        prob_t event_prob(EventClass event_class, event_ind_t event_family, event_ind_t event_row, event_ind_t event_column) const {
            return _vec[this->event_index(event_class, event_family, event_row, event_column)];
        }


        /**
        * \brief Normalise each event family to have sum equal to 1.
        */
        void normaliseEventFamilies() {
            for (event_ind_t i = 2; i < _edges.size(); ++i) {
                prob_t prob_sum = _laplace[i-1] * (_edges[i] - _edges[i-1]);
                for (event_ind_t j = _edges[i-1]; j < _edges[i]; ++j) {
                    prob_sum += _vec[j];
                }

                if (prob_sum) {
                    for (event_ind_t j = _edges[i-1]; j < _edges[i]; ++j) {
                        _vec[j] = (_vec[j] + _laplace[i-1]) / prob_sum;
                    }
                }
            }
        }
        
        void normaliseEventClass(EventClass event_class) {
            event_ind_t event_family = _event_classes[event_class];

            prob_t prob_sum = _laplace[event_family] * (_edges[event_family + 1] - _edges[event_family]);
            for (event_ind_t j = _edges[event_family]; j < _edges[event_family + 1]; ++j) {
                prob_sum += _vec[j];
            }

            if (prob_sum) {
                for (event_ind_t j = _edges[event_family]; j < _edges[event_family + 1]; ++j) {
                    _vec[j] = (_vec[j] + _laplace[event_family]) / prob_sum;
                }
            }
        }


        /**
         *
         */
        event_ind_t max_VJ_ins_len() const {
            return this->eventFamilySize(VJ_VAR_JOI_INS_LEN, 0) - 1;
        }
        event_ind_t max_VD_ins_len() const {
            return this->eventFamilySize(VDJ_VAR_DIV_INS_LEN, 0) - 1;
        }
        event_ind_t max_DJ_ins_len() const {
            return this->eventFamilySize(VDJ_DIV_JOI_INS_LEN, 0) - 1;
        }


        /**
         *
         */
        inline seq_len_t D_min_len(event_ind_t d_index) const { return _d_genes_min_len[d_index - 1]; }


        /**
         * \brief Fill the vector with the given value.
         */
        void fill(prob_t val = 0) {
            _vec[0] = 0;
            for (size_t i = 1; i < _vec.size(); ++i) {
                _vec[i] = val;
            }
        }


        void familyFill(event_ind_t event_family, prob_t val = 0) {
            for (event_ind_t i = _edges[_event_classes[event_family]]; i < _edges[_event_classes[event_family + 1]]; ++i) {
                _vec[i] = val;
            }
        }


        size_t families() const { return _event_classes.size(); }


        size_t size() const { return _vec.size(); }


        seq_len_t n_columns(EventClass event_class, event_ind_t event_family = 0) const {
            return _event_family_col_numbers[_event_classes[event_class] + event_family];
        }

    private:

        vector<prob_t> _vec;
        vector<event_ind_t> _edges;  /** Vector with starting indices for each event family. */
        vector<event_ind_t> _event_classes;  /** Vector of indices of event classes in _edges. */
        vector<seq_len_t> _event_family_col_numbers;  /** Vector of the number of rows of each event family. */
        vector<prob_t> _laplace;
        vector<seq_len_t> _d_genes_min_len;
        Recombination _recomb;


        /**
        * \brief Private default constructor.
        */
        ModelParameterVector() {}

    };
}

#endif