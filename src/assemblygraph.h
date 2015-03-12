/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vadim dot nazarov at mailbox dot com>
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

#ifndef _ASSEMBLYGRAPH_H
#define _ASSEMBLYGRAPH_H


#include "iterator"

#include "matrixchain.h"


namespace ymir {

    struct AssemblyScenarioMatrix;
    class AssemblyGraph;


    /**
    * \struct AssemblyScenarioMatrix
    *
    * \brief Struct for repersenting all assembling scenarios of an immune receptor, i.e. all paths in AssemblyGraph.
    */
    struct AssemblyScenarioMatrix {
        Matrix<prob_t, Dynamic, Dynamic> event_probs;
        Matrix<eventind_t, Dynamic, Dynamic> event_inds;
        vector<numeric> scenario_probs;
        bool normalised;
    };


    struct AssemblyScenarioIterator : public forward_iterator_tag {

    public:

        AssemblyScenarioIterator() {}


        ~AssemblyScenarioIterator() {}


    protected:

        // vector of points with path intersections
        // current selection of points

    };


    /**
    * \class AssemblyGraph
    *
    * \brief Basic class for representing all possible assembling scenarios.
    */
    class AssemblyGraph : protected ProbMatrixChain {

    public:

        /**
        *
        */
        AssemblyGraph() {
            this->_eventind_mat = new EventIndMatrixChain();
        }


        /**
        *
        */
        AssemblyGraph(bool fill_event_indices = true) {
            if (fill_event_indices) {
                this->_eventind_mat = new EventIndMatrixChain();
            } else {
                this->_eventind_mat = nullptr;
            }
        }


        AssemblyGraph(const AssemblyGraph& other) {
            _chain = other._chain;

            _eventind_mat = nullptr;
            if (other._eventind_mat) {
                _eventind_mat = new EventIndMatrixChain(*other._eventind_mat);
            }
        }


        /**
        *
        */
        virtual ~AssemblyGraph() {
            if (this->_eventind_mat) {
                delete this->_eventind_mat;
            }
        }


        /**
        * \brief Make and return matrix with all possible assembly scenarios.
        *
        * \param normalise If true than normalise all scenario probabilities by the full probability.
        *
        * \return AssemblyScenarioMatrix with all assembly scenarios.
        */
        // get not a matrix, but ITERATOR to all scenarios. we don't need all scenarios at once anyway
        const AssemblyScenarioMatrix& allScenarios(bool normalise = true) const {

        }


        /**
        * \brief Compute and return full assembling probability of this sequence.
        *
        * \return Full assembling probability.
        */
        numeric fullProbability() const {
            matrix_t res(this->_chain[0]);
            for (matrix_ind i = 1; i < this->_chain.size(); i++) {
                res = res * this->_chain[i];
            }
            return res(0,0);
        }


        void /* vector of F(k,i) */ forwardAlgorithm() const {

        }


        void /* vector of B(k,i) */ backwardAlgorithm(numeric full_prob) const {

        }


        void baum_welch() {
            // implement it in the Model class

            // compute full probability
            // and vector of fk(i)
            // w/ forward algorithm

            // compute vector of bk(i)
            // w/ backward algorithm

            // multiply by normalised full prob every scenario event
            // and update big vector
            // with scenario event probabilities
        }

    protected:

        EventIndMatrixChain *_eventind_mat;

    };
}

#endif