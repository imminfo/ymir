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

#ifndef _REPERTOIRE_H_
#define _REPERTOIRE_H_


#include "vector"
#include "iterator"

#include "aligner.h"
#include "clonotype.h"
#include "parser.h"


namespace ymir {

    #define DEFAULT_REPERTOIRE_RESERVE_SIZE 300000

    // AssemblyGraphRepertoire aggregate(*ModelParameterVector) - aggregate and compute parameters (bayesian) (add them to this vector)

    /*
     [real data] =>
     => (AbstractParser) => [SequenceRepertoire] =>
     => (AbstractAligner + segments tables) => [AlignmentRepertoire] =>
     => (Model::Indexer) => [IndexedRepertoire] =>
     => (Model::Builder) => [GraphRepertoire]
    */


    /**
    * \class AbstractRepertoireView
    */
    template <class T>
    class AbstractRepertoireView {
    public:

        typedef iterator<random_access_iterator_tag, T> repertoire_iterator;


        AbstractRepertoireView(shared_ptr<vector<T>> pvec, const vector<size_t>& shifts) :
                _source(pvec), _shifts(shifts) {}


        AbstractRepertoireView(const AbstractRepertoireView& other) : _source(other._source), _shifts(other._shifts) {}


        virtual ~AbstractRepertoireView() {}


        virtual const T& operator[] (size_t index) const {
            return this->_source->at(this->_shifts[index]);
        }


        size_t size() const { return this->_shifts.size(); }


    protected:

        shared_ptr<vector<T>> _source;
        vector<size_t> _shifts;


        AbstractRepertoireView() {}

    };


    /**
    * \class AbstractRepertoire
    */
    template <class T>
    class AbstractRepertoire : public AbstractRepertoireView<T> {
    public:


        AbstractRepertoire() {
            vector<T> *_vecobj = new vector<T>();
            shared_ptr<vector<T>> p(_vecobj);
            this->_vec.swap(p);
            this->_vec->reserve(DEFAULT_REPERTOIRE_RESERVE_SIZE);
            this->_source = this->_vec;
            this->_shifts.reserve(DEFAULT_REPERTOIRE_RESERVE_SIZE);
        }


        AbstractRepertoire(const vector<T> &vec) : AbstractRepertoire() {

        }


        virtual ~AbstractRepertoire() {}


        void swap(vector<T>& vec) {
            this->_vec->swap(vec);
            this->_shifts.resize(this->_vec->size());
            for (size_t i = 0; i < this->_vec->size(); ++i) {
                this->_shifts[i] = i;
            }
        }


        T& operator[] (size_t index) {
            return (*(this->_vec))[index];
        }


        AbstractRepertoireView<T> operator[] (const vector<size_t>& indices) {
            return AbstractRepertoireView<T>(this->_vec, indices);
        }


        /**
        * \brief Get read-only access to the first N elements.
        *
        * \param size N first elements.
        *
        * \return AbstractRepertoireView to the first N elements.
        */
        AbstractRepertoireView<T> head(size_t size) const {
            vector<size_t> shifts;
            shifts.reserve(size);
            for (size_t i = 0; i < size; ++i) {
                shifts.push_back(i);
            }
            return (*this)[shifts];
        }


        /**
        * \brief Get read-only access to the sub-vector - continous sample of the repertoire.
        *
        * \param start Start index of the sub-vector.
        * \param end End index of the sub-vector.
        *
        * \return AbstractRepertoireView to the sub-vector.
        */
        AbstractRepertoireView<T> slice(size_t start, size_t end) const {
            vector<size_t> shifts;
            shifts.reserve(end - start + 1);
            for (size_t i = start; i < end; ++i) {
                shifts.push_back(i);
            }
            return (*this)[shifts];
        }


        /**
        * \brief Get read-only access to a random sample of the repertoire.
        *
        * \param size Size of the sample.
        * \param seed Seed for generating random numbers.
        * \param replace Should sampling be with replacement or not.
        *
        * \return AbstractRepertoireView to the random sample.
        */
        AbstractRepertoireView<T> sample(size_t size, size_t seed, bool replace = true) const {
            // !!
            // !!
            // !!
        }


    protected:

        shared_ptr<vector<T>> _vec;


        AbstractRepertoire(const AbstractRepertoire<T>& other) {}

    };


//    typedef AbstractRepertoireView<Clonotype> ClonesetView;

    class ClonesetView : public AbstractRepertoireView<Clonotype> {

    public:

        ClonesetView() {

        }


        virtual ~ClonesetView() {

        }


//        ClonesetView codingSequences() const {
//
//        }


//        ClonesetView noncodingSequences() const {
//
//        }

    };


    typedef AbstractRepertoire<Clonotype> Cloneset;

}

#endif