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
using namespace Eigen;


namespace ymir {

    template <typename _Scalar>
    class MatrixChain;


    /**
    * \class EventIndMatrixChain
    *
    * \brief Class for storing chain of matrices of indices of scenario events.
    */
    typedef MatrixChain<eventind_t> EventIndMatrixChain;


    /**
    * \class ProbMatrixChain
    *
    * \brief Class for storing chain of matrices of scenario event probabilities.
    */
    typedef MatrixChain<prob_t> ProbMatrixChain;


    /**
    * \class MatrixChain
    *
    * \brief Basic class for representing chains (lists) of matrices, i.e. graph representation for receptors or,
    * more generally, for probabilistic graphical models.
    *
    */
    template <typename _Scalar>
    class MatrixChain {

    protected:

        /**
        * \typedef matrix_t
        *
        * \brief Type of matrices in the chain.
        */
        typedef Matrix<_Scalar, Dynamic, Dynamic> matrix_t;


    public:

        /**
        * \typedef matrix_ind
        *
        * \brief Index type for accessing matrices in the chain.
        */
        typedef uint8_t matrix_ind;


        /**
        * \brief Basic constructor, adds first 1x1 matrix (starting point in
        * a probabilistic graphical model)
        */
        MatrixChain() {
            this->addMatrix(1, 1);
        }


        MatrixChain(const MatrixChain& other) {
            _chain = other._chain;
        }


        /**
        *
        */
        virtual ~MatrixChain() {}


        /**
        * \brief Add new matrix with specified rows and cols.
        *
        * \param nrow Number of rows in the new matrix.
        * \param ncol Number of cols in the new matrix.
        *
        * \return Index in the chain of the new matrix.
        */
        matrix_ind addMatrix(size_t nrow, size_t ncol) {
            this->_chain.push_back(matrix_t(nrow, ncol));
            return this->_chain.size() - 1;
        }


        /**
        * \brief Access element with specified indices at specified matrix.
        *
        * \param i Row of the element.
        * \param j Column of the element.
        * \param mat_i Element's matrix.
        *
        * \return Element at position (i,j) from matrix mat_i.
        */
        const _Scalar& operator()(size_t i, size_t j, matrix_ind mat_i) const {
            return this->_chain[mat_i](i, j);
        }

        _Scalar& operator()(size_t i, size_t j, matrix_ind mat_i) {
            return this->_chain[mat_i](i, j);
        }


        /**
        * \brief Get the number of matrices in this chain.
        */
        matrix_ind matrices() const {
            return this->_chain.size();
        }


        /**
        * \brief Get the number of rows in a matrix with the given index.
        */
        size_t rows(matrix_ind mat_i) const {
            return this->_chain[mat_i].rows();
        }


        /**
        * \brief Get the number of columns in a matrix with the given index.
        */
        size_t cols(matrix_ind mat_i) const {
            return this->_chain[mat_i].cols();
        }


        /**
        * \brief Iterative product of all matrices in the chain.
        */
        _Scalar chainProd(bool *err = nullptr) const {
            matrix_t res(this->_chain[0]);

            if (err) {
                *err = false;
            }

            for (matrix_ind i = 1; i < this->_chain.size(); i++) {
                if (res.cols() != this->_chain[i].rows()) {
                    if (err) {
                        *err = true;
                    }
                    break;
                }
                res = res * this->_chain[i];
            }
            return res(0,0);
        }


    protected:

        vector<matrix_t> _chain;
    };
}