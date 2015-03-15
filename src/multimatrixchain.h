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

#ifndef _MULTIMATRIXCHAIN_H
#define _MULTIMATRIXCHAIN_H


#include <vector>

#include "types.h"


using namespace std;
using namespace Eigen;


namespace ymir {

    template <typename _Scalar>
    class MultiMatrixChain;


    template <typename _Scalar>
    class MMCSlice;


    /**
    * \class EventIndMatrixChain
    *
    * \brief Class for storing chain of matrices of indices of scenario events.
    */
    typedef MultiMatrixChain<eventind_t> EventIndMMC;

    typedef MMCSlice<eventind_t> EventIndMMCSlice;


    /**
    * \class ProbMatrixChain
    *
    * \brief Class for storing chain of matrices of scenario event probabilities.
    */
    typedef MultiMatrixChain<prob_t> ProbMMC;

    typedef MMCSlice<prob_t> ProbMMCSlice;


    /**
    * \class MultiMatrixChain
    *
    * \brief Class for storing lists of matrices, where one node in the list (called "chain") could
    * contain more than one matrix.
    */
    template <typename _Scalar>
    class MultiMatrixChain {

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
        * \typedef matrix_t
        *
        * \brief Type of matrices in the chain.
        */
        typedef Matrix<_Scalar, Dynamic, Dynamic> matrix_t;


        /**
        * \typedef dim_t
        *
        * \brief Type of dimensions of matrices (rows and columns).
        */
        typedef size_t dim_t;

    protected:

        /**
        * \struct Node
        *
        * \brief Node in the chain. Stores one or more matrices with equal size.
        */
        struct Node {
        public:

            Node(matrix_ind_t n, dim_t rows, dim_t cols) : _n(n) {
                _vec = new matrix_t[n];
                for (matrix_ind_t i = 0; i < _n; ++i) {
                    _vec[i].resize(rows, cols);
                    _vec[i].fill(0);
                }
                // _vec = new matrix_t*[n];
                // _vec[i] = new matrix_t(rows, cols);
            }


            ~Node() { delete [] _vec; }


            ///@{
            const matrix_t& operator[](matrix_ind_t mat_i) const { return _vec[mat_i]; }
            matrix_t& operator[](matrix_ind_t mat_i) { return _vec[mat_i]; }
            ///@}


            matrix_ind_t size() const { return _n; }

//            matrix_t *ptr(matrix_ind_t mat_i) const { return _vec + mat_i; }

        protected:

            matrix_t *_vec;
            matrix_ind_t _n;


            Node() {}

        };

    public:

        MultiMatrixChain() {

        }


        virtual ~MultiMatrixChain() {
        }


        /**
        * \brief Access to matrices in nodes.
        */
        const matrix_t& matrix(node_ind_t node_i, matrix_ind_t mat_i) const {
            return _chain[node_i][mat_i];
        }


        node_ind_t chainSize() const { return _chain.size(); }


        matrix_ind_t nodeSize(node_ind_t node_i) const { return _chain[node_i].size(); }


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
        const _Scalar& operator()(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) const {
            return _chain[node_i][mat_i](row, col);
        }
        _Scalar& operator()(node_ind_t node_i, matrix_ind_t mat_i, dim_t row, dim_t col) {
            return _chain[node_i][mat_i](row, col);
        }
        ///@}


        /**
        * \brief Add new node with pattern matrix.
        */
        node_ind_t addNode(node_ind_t node_i, matrix_ind_t n_matrices, dim_t rows, dim_t cols) {
            _chain.push_back(Node(n_matrices, rows, cols));
            return _chain.size() - 1;
        }


//        MMCSlice<_Scalar> slice() const {
//            // ???
//            // ???
//            // ???
//        }


    protected:

        vector<Node> _chain;

    };


    /**
    * \class MMCSlice
    *
    * \brief Plain list of matrices, i.e., without multiple matrices at each node.
    */
    template <typename _Scalar>
    class MMCSlice {
    friend class MultiMatrixChain<_Scalar>; // yes-yes, I know friends are bad.

    public:

        MMCSlice(typename MultiMatrixChain<_Scalar>::node_ind_t n) : _i(0), _n(n) {
            this->_vec = new MultiMatrixChain<_Scalar>::matrix_t*[_n];
        }


        ~MMCSlice() { delete [] this->_vec; }


        const typename MultiMatrixChain<_Scalar>::matrix_t& operator[](typename MultiMatrixChain<_Scalar>::node_ind_t mat_i) const {
            return *(this->_vec[mat_i]);
        }

    protected:

        typename MultiMatrixChain<_Scalar>::matrix_t** _vec;
        typename MultiMatrixChain<_Scalar>::node_ind_t _i, _n;


        MMCSlice<_Scalar>& addMatrix(typename MultiMatrixChain<_Scalar>::matrix_t *pmat) {
            this->_vec[_i] = pmat;
            ++_i;
            return *this;
        }
    };
}

#endif