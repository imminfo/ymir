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


#include <vector>

#include "types.h"


using namespace std;


namespace ymir {

    /**
    * \class ModelParameterVector
    *
    * \brief Class for storing parameters of assembling statistical model. Note:
    * event with index 0 (zero) is "null" event and always has zero probability.
    */
    class ModelParameterVector {
    public:

        /**
        * \brief Construct the vector from the given vector of event probabilities and lengths of each event family (V deletions, insertions length, etc.).
        *
        * \param param_vec Vector with event probabilities.
        * \param lens_vec Vector with lengths of each event family.
        * \param laplace_vec Vector with laplace correction value for each event family.
        * \param do_normalise Boolean, do normalisation of event families after initialising?
        */
        ModelParameterVector(const vector<prob_t>& param_vec, const vector<eventind_t>& lens_vec,
                const vector<prob_t>& laplace_vec = vector<prob_t>(), bool do_normalise = true)
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

    private:

        vector<prob_t> _vec;
        vector<eventind_t> _edges;
        vector<prob_t> _laplace;


        /**
        * \brief Private default constructor.
        */
        ModelParameterVector() {}

    };
}