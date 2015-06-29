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

#ifndef _MULTIMATRIXCHAIN_H_
#define _MULTIMATRIXCHAIN_H_


#include <vector>

#include "types.h"


#define YDEBUG

using namespace std;


namespace ymir {

    template <typename _Scalar>
    class MultiMatrixChain;


    /**
    * \class EventIndMatrixChain
    *
    * \brief Class for storing chain of matrices of indices of scenario events.
    */
    typedef MultiMatrixChain<eventind_t> EventIndMMC;


    /**
    * \class ProbMatrixChain
    *
    * \brief Class for storing chain of matrices of scenario event probabilities.
    */
    typedef MultiMatrixChain<prob_t> ProbMMC;


    /**
    * \class MultiMatrixChain
    *
    * \brief Class for storing lists of matrices, where one node in the list (called "chain") could
    * contain more than one matrix.
    */
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
//        typedef Matrix<_Scalar, Dynamic, Dynamic> matrix_t;
        typedef Matrix<_Scalar, dim_t> matrix_t;

    protected:

        /**
        * \struct Node
        *
        * \brief Node in the chain. Stores one or more matrices with equal size.
        */
        struct Node {
        public:

            Node() : _n(0), _start_index(0), _rows(0), _cols(0) {}


            Node(size_t start_index, matrix_ind_t n, dim_t rows, dim_t cols) {
                this->init(start_index, n, rows, cols);
            }


            Node(const Node& other) {
                _n = other._n;
                _rows = other._rows;
                _cols = other._cols;
                _start_index = other._start_index;
            }


            ~Node() { }


            Node& operator= (const Node& other) {
                _n = other._n;
                _rows = other._rows;
                _cols = other._cols;
                _start_index = other._start_index;
                return *this;
            }


            void init(size_t start_index, matrix_ind_t n, dim_t rows, dim_t cols) {
                _n = n;
                _rows = rows;
                _cols = cols;
                _start_index = start_index;
            }


            size_t operator() (matrix_ind_t mat, dim_t row, dim_t column) const {
#ifdef YDEBUG
                if (!(row >= 0 && row < _rows && column >= 0 && column < _cols)) { throw(std::runtime_error("Rows / columns number check failed!")); }
#endif
                return _start_index + mat * (_rows * _cols) + row * _cols + column;
            }


            matrix_ind_t size() const { return _n; }

            dim_t rows() const { return _rows; }

            dim_t cols() const { return _cols; }

            size_t start() const { return _start_index; }


        protected:

            size_t _start_index;
            matrix_ind_t _n;
            dim_t _rows, _cols;

        };

    public:

        MultiMatrixChain() {
            _chain.clear();
            _values.clear();
            _values.reserve(5000);
        }


        MultiMatrixChain(const MultiMatrixChain &other) {
            _chain = other._chain;
            _values = other._values;
        }


        MultiMatrixChain(MultiMatrixChain &&other) {
            _chain.swap(other._chain);
            _values.swap(other._values);
        }


        virtual ~MultiMatrixChain() { }


        MultiMatrixChain& operator= (const MultiMatrixChain &other) {
            _chain = other._chain;
            _values = other._values;
            return *this;
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
            return _values[_chain[node_i](mat_i, row, col)];
        }

        const _Scalar& operator()(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
            return _values[_chain[node_i](mat_i, row, col)];
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
            initNode(_chain.size() - 1, n_matrices, rows, cols);
            return _chain.size() - 1;
        }
        ///@}


        void initNode(node_ind_t node_i, matrix_ind_t n_matrices, dim_t rows, dim_t cols) {
            _chain[node_i].init(_values.size(), n_matrices, rows, cols);
            _values.resize(_values.size() + rows * cols * n_matrices, 0);
        }


        void swap(MultiMatrixChain<_Scalar> &other) {
            _chain.swap(other._chain);
            _values.swap(other._values);
        }


        void finish() {
            _chain.shrink_to_fit();
            _values.shrink_to_fit();
//            cout << _values.size() << endl;
//            cout << _values.capacity() << endl;
        }


        void fill(node_ind_t node, matrix_ind_t mat, _Scalar val = 0) {
            for (dim_t r = 0; r < _chain[node].rows(); ++r) {
                for (dim_t c = 0; c < _chain[node].cols(); ++c) {
                    _values[_chain[node](mat, r, c)] = val;
                }
            }
        }


        matrix_t matrix(node_ind_t node, matrix_ind_t mat) const {
            matrix_t res(_chain[node].rows(), _chain[node].cols(), 0);
            for (dim_t r = 0; r < _chain[node].rows(); ++r) {
                for (dim_t c = 0; c < _chain[node].cols(); ++c) {
                    res(r, c) = (*this)(node, mat, r, c);
                }
            }
            return res;
        }


        _Scalar operator[](size_t index) const { return _values[index]; }

        size_t values_size() const { return _values.size(); }

    protected:

        vector<Node> _chain;
        vector<_Scalar> _values;

    };

}

#endif