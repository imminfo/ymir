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

/*
 *  MULTI MATRIX CHAIN W/ EIGEN
 */

#ifndef _MULTIMATRIXCHAIN_H_
#define _MULTIMATRIXCHAIN_H_


#include <vector>

#include "types.h"
#include "Eigen/Core"
#include "Eigen/SparseCore"


//#define YDEBUG

namespace ymir {

    template <typename _Scalar>
    class MultiMatrixChain;


    /**
    * \class EventIndMatrixChain
    *
    * \brief Class for storing chain of matrices of indices of scenario events.
    */
    typedef MultiMatrixChain<event_ind_t> EventIndMMC;

    typedef unique_ptr<EventIndMMC> pEventIndMMC;


    /**
    * \class ProbMatrixChain
    *
    * \brief Class for storing chain of matrices of scenario event probabilities.
    */
    typedef MultiMatrixChain<prob_t> ProbMMC;

    typedef unique_ptr<ProbMMC> pProbMMC;


    /**
     * \class NumErrorsMMC
     *
     * \brief Class for storing a number of errors for each scenario event.
     */
    typedef seq_len_t error_num_t;
    typedef MultiMatrixChain<error_num_t> ErrMMC;

    typedef unique_ptr<ErrMMC> pErrMMC;


    typedef MultiMatrixChain<codon_hash> CodonMMC;


    template <typename _Scalar>
    class MultiMatrixChain {

        friend class MAAGForwardBackwardAlgorithm;  // ):

    public:

        /**
        * \typedef node_ind_t
        *
        * \brief Node index type.
        */
        typedef uint8_t node_ind_t;


        /**
        * \typedef matrix_ind_t
        *
        * \brief Matrix index type.
        */
        typedef uint8_t matrix_ind_t;


        /**
        * \typedef dim_t
        *
        * \brief Type of dimensions of matrices (rows and columns).
        */
//        typedef size_t dim_t;
        typedef seq_len_t dim_t;


        /**
        * \typedef matrix_t
        *
        * \brief Type of matrices in the chain.
        */
        typedef Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> matrix_t;

    protected:

        /**
         * \struct Node
         *
         * \brief Node in the chain. Stores one or more matrices with equal size.
         */
        struct Node {
        public:

            Node() : _n(0), _rows(0), _cols(0)
            { }


            void init(size_t start_index, matrix_ind_t n, dim_t rows, dim_t cols) {
                _n = n;
                _rows = rows;
                _cols = cols;
                for (int i = 0; i < _n; ++i) {
                    _matrices.push_back(matrix_t::Zero(_rows, _cols));
                }
            }


            _Scalar operator() (matrix_ind_t mat, dim_t row, dim_t column) const {
#ifndef DNDEBUG
                if (!(mat >= 0 && mat < _n)) { throw(std::runtime_error("Matrix index check failed!")); }

                if (!(row >= 0 && row < _rows && column >= 0 && column < _cols)) { throw(std::runtime_error("Rows / columns number check failed!")); }
#endif
                return _matrices[mat](row, column);
            }


            _Scalar& operator() (matrix_ind_t mat, dim_t row, dim_t column) {
#ifndef DNDEBUG
                if (!(mat >= 0 && mat < _n)) { throw(std::runtime_error("Matrix index check failed!")); }

                if (!(row >= 0 && row < _rows && column >= 0 && column < _cols)) { throw(std::runtime_error("Rows / columns number check failed!")); }
#endif
                return _matrices[mat](row, column);
            }


            matrix_ind_t size() const { return _n; }

            dim_t rows() const { return _rows; }

            dim_t cols() const { return _cols; }

//            size_t start() const { return _start_index; }

            size_t n_values() const { return _n * _rows * _cols; }

            const matrix_t& matrix(matrix_ind_t mat) const {
#ifndef DNDEBUG
                if (!(mat >= 0 && mat < _n)) { throw(std::runtime_error("Matrix index check failed!")); }
#endif
                return _matrices[mat];
            }


        protected:

            matrix_ind_t _n;
            std::vector<matrix_t> _matrices;
            dim_t _rows, _cols;

        };

    public:

        MultiMatrixChain() {
        }


//        MultiMatrixChain(MultiMatrixChain &) = default;
//
//
//        MultiMatrixChain(MultiMatrixChain &&) = default;
//
//
//        MultiMatrixChain& operator=(MultiMatrixChain &) = default;
//
//
//        MultiMatrixChain& operator=(MultiMatrixChain &&) = default;


        virtual ~MultiMatrixChain()
        {
        }


        void resize(node_ind_t n_nodes) {
            _chain.resize(n_nodes);
        }


        node_ind_t chainSize() const { return _chain.size(); }


        ///@{
        matrix_ind_t nodeSize(node_ind_t node_i) const { return _chain[node_i].size(); }

        dim_t nodeRows(node_ind_t node_i) const { return _chain[node_i].rows(); }

        dim_t nodeColumns(node_ind_t node_i) const { return _chain[node_i].cols(); }
        ///@}


        /**
        * \brief Access element with specified indices at specified matrix.
        *
        * \param node_i Node's index with this matrix.
        * \param mat_i Element's matrix.
        * \param row Row of the element.
        * \param col Column of the element.
        *
        * \return Element at position (i,j) from matrix mat_i from node node_i.
        */
        ///@{
        _Scalar& operator()(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) {
#ifndef DNDEBUG
            if (node_i >= _chain.size()) { throw(std::runtime_error("Number of the Node is out of bounds."));}
#endif
            return _chain[node_i](mat_i, row, col);
        }

//        const _Scalar& operator()(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
        _Scalar operator()(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
#ifndef DNDEBUG
            if (node_i >= _chain.size()) { throw (std::runtime_error("Number of the Node is out of bounds.")); }
#endif
            return _chain[node_i](mat_i, row, col);
        }

        _Scalar& at(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) {
#ifndef DNDEBUG
            if (node_i >= _chain.size()) { throw(std::runtime_error("Number of the Node is out of bounds."));}
#endif
            return _chain[node_i](mat_i, row, col);
        }

//        const _Scalar& operator()(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
        _Scalar at(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
#ifndef DNDEBUG
            if (node_i >= _chain.size()) { throw(std::runtime_error("Number of the Node is out of bounds."));}
#endif
            return _chain[node_i](mat_i, row, col);
        }
        ///@}


        /**
        * \brief Add new node with pattern matrix.
        */
        ///@{
        node_ind_t addNode() {
            _chain.push_back(Node());
            return _chain.size() - 1;
        }

        node_ind_t addNode(matrix_ind_t n_matrices, dim_t rows, dim_t cols) {
            _chain.push_back(Node());
            this->initNode(_chain.size() - 1, n_matrices, rows, cols);
            return _chain.size() - 1;
        }
        ///@}


        void initNode(node_ind_t node_i, matrix_ind_t n_matrices, dim_t rows, dim_t cols) {
            size_t prev_node_size = _chain[node_i].n_values();
            _chain[node_i].init(_values.size(), n_matrices, rows, cols);
//            _values.resize(_values.size() + _chain[node_i].n_values() - prev_node_size, 0);
        }


        void swap(MultiMatrixChain &other) {
            _chain.swap(other._chain);
        }


        void finish() {
            _chain.shrink_to_fit();
//            cout << _values.size() << endl;
//            cout << _values.capacity() << endl;
        }


        void fill(node_ind_t node, matrix_ind_t mat, _Scalar val = 0) {
            for (dim_t r = 0; r < _chain[node].rows(); ++r) {
                for (dim_t c = 0; c < _chain[node].cols(); ++c) {
                    _chain[node](mat, r, c) = val;
                }
            }
        }


        const matrix_t& matrix(node_ind_t node, matrix_ind_t mat) const {
            return _chain[node].matrix(mat);
        }


//        _Scalar operator[](size_t index) const { return _values[index]; }
//        size_t values_size() const { return _values.size(); }


    protected:

        std::vector<Node> _chain;
        std::vector<_Scalar> _values;

    };

}

#endif