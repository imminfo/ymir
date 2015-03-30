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
        * \brief Construct the vector from the given vector of event probabilities and lengths of each event family (V deletions, insertions length, etc.).
        * for the VJ-recombination model.
        *
        * \param param_vec Vector with event probabilities.
        * \param lens_vec Vector with lengths of each event family.
        * \param laplace_vec Vector with laplace correction value for each event family.
        * \param do_normalise Boolean, do normalisation of event families after initialising?
        */
        ModelParameterVector(const vector<prob_t>& param_vec,
                const vector<eventind_t>& lens_vec,
                const vector<prob_t>& laplace_vec = vector<prob_t>(),
                bool do_normalise = true)
        {
            _vec = vector<prob_t>();
            _vec.reserve(param_vec.size() + 1);
            _vec.push_back(0);
            for (eventind_t i = 0; i < param_vec.size(); ++i) {
                _vec.push_back(param_vec[i]);
            }

            _edges.reserve(lens_vec.size() + 2);

            _edges.push_back(0);
            _edges.push_back(1);
            for (eventind_t i = 0; i < lens_vec.size(); ++i) {
                _edges.push_back(_edges[i+1] + lens_vec[i]);
            }
            _edges.push_back(_vec.size());

            _laplace = vector<prob_t>();
            if (laplace_vec.size()) {
                _laplace.reserve(laplace_vec.size());
                _laplace.push_back(0);
                for (eventind_t i = 0; i < _edges.size() - 1; ++i) {
                    _laplace.push_back(laplace_vec[i]);
                }
            } else {
                // Laplace correction equal to zero if vector is not supplied.
                _laplace.resize(_edges.size() - 1, 0);
            }

            if (do_normalise) {
                this->normaliseEventFamilies();
            }
        }


        /**
        * \brief Construct the vector from the given vector of event probabilities and lengths of each event family (V deletions, insertions length, etc.).
        * for the VDJ-recombination model.
        *
        * \param param_vec Vector with event probabilities.
        * \param lens_vec Vector with lengths of each event family.
        * \param d_gene_max_dels Vector with number of maximal number of D'3 deletions for each Diversity segment.
        * \param laplace_vec Vector with laplace correction value for each event family.
        * \param do_normalise Boolean, do normalisation of event families after initialising?
        */
        ModelParameterVector(const vector<prob_t>& param_vec,
                const vector<eventind_t>& lens_vec,
                vector<seq_len_t> d_gene_max_dels,
                const vector<prob_t>& laplace_vec = vector<prob_t>(),
                bool do_normalise = true) : ModelParameterVector(param_vec, lens_vec, laplace_vec, do_normalise)
        {

            _d_gene_num = d_gene_max_dels.size();
            _d_gene_max_dels = d_gene_max_dels;
        }


        //============= VECTOR INDICES ACCESS =============//


        /**
        * \brief Get a probability of an event with the given event family's index and event's local index.
        *
        * \param event_family Event family's index.
        * \param event_index Local index of the event.
        *
        * \return Probability of the event.
        */
        const prob_t& getEventProbability(eventind_t event_family, eventind_t event_index) const {
            return _vec[_edges[event_family] + event_index];
        }


        eventind_t eventFamilySize(eventind_t event_family) const {
            return _edges[event_family + 1] - _edges[event_family];
        }


        /**
        * \brief Get a probability of an event with the given global index.
        *
        * \param event_index Global index of the event.
        *
        * \return Probability of the event.
        */
        const prob_t& getEventProbability(eventind_t event_index) const {
            return _vec[event_index];
        }
        const prob_t& operator[] (eventind_t event_index) const {
            return _vec[event_index];
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


        vector<prob_t>::const_iterator get_iterator(eventind_t i) const {
            return _vec.begin() + i;
        }


        //============= EVENT ACCESS =============//


        inline eventind_t index_V_gene(segindex_t v_index) const {
            return _edges[1] + v_index - 1;
        }
        prob_t prob_V_gene(segindex_t v_index) const { return _vec[index_V_gene(v_index)]; }  // Hmmm...


        inline eventind_t index_J_gene(segindex_t j_index) const { return _edges[2] + j_index - 1; }
        prob_t prob_J_gene(segindex_t j_index) const { return _vec[index_J_gene(j_index)]; }


        inline eventind_t index_JD_genes(segindex_t j_index, segindex_t d_index) const { return _edges[2] + (j_index - 1) * _d_gene_num + (d_index - 1); }
        prob_t prob_JD_genes(segindex_t j_index, segindex_t d_index) const { return _vec[index_JD_genes(j_index, d_index)]; }


        inline eventind_t index_V_del(segindex_t v_index, seq_len_t del_num) const { return _edges[3 + (v_index - 1)] + del_num; }
        prob_t prob_V_del(segindex_t v_index, seq_len_t del_num) const { return _vec[index_V_del(v_index, del_num)]; }


        inline eventind_t index_J_del(segindex_t j_index, seq_len_t del_num) const { return _edges[3 + (_edges[2] - 1) + (j_index - 1)] + del_num; }
        prob_t prob_J_del(segindex_t j_index, seq_len_t del_num) const { return _vec[index_J_del(j_index, del_num)]; }


        inline eventind_t index_D_del(segindex_t d_index, seq_len_t d5_del_num, seq_len_t d3_del_num) const {
            return _edges[3 + (_edges[2] - 1) + 2 + (d_index - 1)] + d5_del_num*_d_gene_max_dels[d_index - 1] + d3_del_num;
        }
        prob_t prob_D_del(segindex_t d_index, seq_len_t d5_del_num, seq_len_t d3_del_num) const {
            return _vec[index_D_del(d_index, d5_del_num, d3_del_num)];
        }


        inline eventind_t index_VJ_ins_len(seq_len_t ins_len) const {
            return _edges[_edges.size() - 7] + ins_len; // - 1 - 4 - 2
        }
        prob_t prob_VJ_ins_len(seq_len_t ins_len) const {
            return _vec[index_VJ_ins_len(ins_len)];
        }
        eventind_t max_VJ_ins_len() const {
            return this->eventFamilySize(_edges.size() - 7) - 1;
        }


        inline eventind_t index_VD_ins_len(seq_len_t ins_len) const {
            return _edges[_edges.size() - 12] + ins_len; // - 1 - 8 - 3
        }
        prob_t prob_VD_ins_len(seq_len_t ins_len) const {
            return _vec[index_VD_ins_len(ins_len)];
        }
        eventind_t max_VD_ins_len() const {
            return this->eventFamilySize(_edges.size() - 12) - 1;
        }


        inline eventind_t index_DJ_ins_len(seq_len_t ins_len) const {
            return _edges[_edges.size() - 11] + ins_len; // -1 - 8 - 2
        }
        prob_t prob_DJ_ins_len(seq_len_t ins_len) const {
            return _vec[index_DJ_ins_len(ins_len)];
        }
        eventind_t max_DJ_ins_len() const {
            return this->eventFamilySize(_edges.size() - 11) - 1;
        }


        inline eventind_t index_VJ_ins_nuc() const {
            return _edges[_edges.size() - 6]; // -1 - 4 - 1
        }


        inline eventind_t index_VD_ins_nuc() const {
            return _edges[_edges.size()  - 10]; // -1 - 8 - 1
        }


        inline eventind_t index_DJ_ins_nuc() const {
            return _edges[_edges.size() - 6]; // -1 - 4 - 1
        }

    private:

        vector<prob_t> _vec;
        vector<eventind_t> _edges;  /** Vector with starting indices for each event family. */
        vector<prob_t> _laplace;
        segindex_t _d_gene_num;
        vector<seq_len_t> _d_gene_max_dels;


        /**
        * \brief Private default constructor.
        */
        ModelParameterVector() {}

    };
}

#endif