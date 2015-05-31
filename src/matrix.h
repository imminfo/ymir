//
// Created by Vadim N. on 27/05/2015.
//

#ifndef YMIR_MATRIX_H
#define YMIR_MATRIX_H


//#include "types.h"
#include <stdexcept>


#define YDEBUG

namespace ymir {


    template <typename _Scalar, typename _Dim>
    class Matrix {

    public:

        // def
        Matrix() {
            _rows = 0;
            _cols = 0;
            _data = nullptr;
        }

        Matrix(_Dim rows, _Dim columns, _Scalar val = 0) : _rows(rows), _cols(columns) {
            _data = new _Scalar[_rows * _cols];
            this->fill(val);
        }

        // copy
        Matrix(const Matrix &other) {
            _rows = other._rows;
            _cols = other._cols;
            if (_data) { delete [] _data; }
            if (_rows * _cols != 0) {
                _data = new _Scalar[_rows * _cols];
                for (_Dim r = 0; r < _rows; ++r) {
                    for (_Dim c = 0; c < _cols; ++c) {
                        _data[r * _cols + c] = other._data[r * _cols + c];
                    }
                }
            } else {
                _data = nullptr;
            }
        }

        // dest
        virtual ~Matrix() {
            if (_data) { delete [] _data; }
        }


        void fill(_Scalar val = 0) {
            for (_Dim r = 0; r < _rows; ++r) {
                for (_Dim c = 0; c < _cols; ++c) {
                    _data[r * _cols + c] = val;
                }
            }
        }

        // eq
        Matrix& operator=(const Matrix &other) {
            _rows = other._rows;
            _cols = other._cols;
            if (_data) { delete [] _data; }
            if (_rows * _cols != 0) {
                _data = new _Scalar[_rows * _cols];
                for (_Dim r = 0; r < _rows; ++r) {
                    for (_Dim c = 0; c < _cols; ++c) {
                        _data[r * _cols + c] = other._data[r * _cols + c];
                    }
                }
            } else {
                _data = nullptr;
            }
            return *this;
        }

        // resize
        void resize(_Dim rows, _Dim columns) {
            _rows = rows;
            _cols = columns;
            if (_data) { delete [] _data; }
            _data = new _Scalar[_rows * _cols];
        }

        // rows
        _Dim rows() const { return _rows; }

        // cols
        _Dim cols() const { return _cols; }


// operator() with check
        _Scalar& operator()(_Dim row, _Dim col) {
#ifdef YDEBUG
            if (!(row >= 0 && row < _rows && col >= 0 && col < _cols)) { throw(std::runtime_error("Rows / columns number check failed!")); }
#endif
            return _data[row * _cols + col];
        }

        const _Scalar& operator()(_Dim row, _Dim col) const {
#ifdef YDEBUG
            if (!(row >= 0 && row < _rows && col >= 0 && col < _cols)) { throw(std::runtime_error("Rows / columns number check failed!")); }
#endif
            return _data[row * _cols + col];
        }

        // operator ==
        bool operator==(const Matrix &other) const {
            if (_rows != other._rows || _cols != other._cols) { return false; }

            for (_Dim r = 0; r < _rows; ++r) {
                for (_Dim c = 0; c < _cols; ++c) {
                    if (_data[r * _cols + c] != other._data[r * _cols + c]) {
                        return false;
                    }
                }
            }

            return true;
        }

        // operator*
        Matrix operator*(const Matrix &other) const {
#ifdef YDEBUG
            if (_cols != other._rows) { throw(std::runtime_error("Multiplication of matrices with wrong dimensions!")); }
#endif

            Matrix res(_rows, other._cols, 0);
            for (_Dim i = 0; i < _rows; ++i) {
                for (_Dim j = 0; j < other._cols; ++j) {
                    for (_Dim k = 0; k < _cols; ++k) {
                        res(i, j) += (*this)(i, k) * other(k, j);
                    }
                }
            }
            return res;
        }


        Matrix operator*(const _Scalar &val) const {
            Matrix res(_rows, _cols);
            for (_Dim i = 0; i < _rows; ++i) {
                for (_Dim j = 0; j < _cols; ++j) {
                    res(i, j) = res(i, j) * val;
                }
            }
            return res;
        }


    protected:

        _Scalar *_data;
        _Dim _rows, _cols;

    };


    template <typename _Scalar, typename _Dim>
    inline Matrix<_Scalar, _Dim> operator*(_Scalar val, const Matrix<_Scalar, _Dim> &rhs) {
        Matrix<_Scalar, _Dim> res(rhs.rows(), rhs.cols());
        for (_Dim i = 0; i < rhs.rows(); ++i) {
            for (_Dim j = 0; j < rhs.cols(); ++j) {
                res(i, j) = res(i, j) * val;
            }
        }
        return res;
    }

}

#endif //YMIR_MATRIX_H
