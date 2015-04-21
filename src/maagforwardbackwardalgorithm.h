//
// Created by Vadim N. on 20/04/2015.
//

#ifndef YMIR_MAAGFORWARDBACKWARDALGORITHM_H
#define YMIR_MAAGFORWARDBACKWARDALGORITHM_H

#include "maag.h"

namespace ymir {


    class MAAGForwardBackwardAlgorithm : protected MAAG {
    public:


        MAAGForwardBackwardAlgorithm(const MAAG &maag) {
            init(maag);
        }


        virtual ~MAAGForwardBackwardAlgorithm() { }


        bool init(const MAAG &maag) {
            _status = false;
            _node_i = 0;
            _mat_i = 0;
            _row = 0;
            _column = 0;
            if (maag && maag->_events) {
                _chain = maag._chain;
                _events = maag._events;
                _status = true;
                _recomb = maag.recombination();
            }
            return _status;
        }


        event_pair_t next() {
            if (_status && _node_i < _events->chainSize()) {

            }

            _status = false;
            retrun event_pair_t(0, 0);
        }


        bool status() const { return _status; }

    protected:

        bool _status;
        RECOMBINATION _recomb;
        size_t _node_i, _mat_i, _row, _column;
        // _memo;


        prob_t forward(node_ind_t node_i,
                       dim_t row_i,
                       dim_t column_i,
                       const vector<matrix_ind_t> &matrix_indices) {
            // find previous nodes, compute forward probabilities for them, sum them up and return the sum
            prob_t res = 1;
            node_ind_t prev_node_i = node_i - 1;
            if (prev_node_i > 0 && node_i < _chain.size() - 1) {
                for (dim_t next_row_i = 0; next_row_i < nodeRows(prev_node_i); ++next_row_i) {
                    res += (*this)(node_i, matrix_indices[node_i], row_i, column_i) * forward(prev_node_i, matrix_indices[prev_node_i], next_row_i, row_i);
                }
            } else if (prev_node == 0) {
                if (_recomb == VJ_RECOMB) {
                    res = _chain[0](matrix_indices[0], matrix_indices[VJ_CHAIN_SIZE - 1]);
                } else {
                    res = _chain[0][matrix_indices[0]](0, 0);
                }
            } else {
                if (_recomb == VJ_RECOMB) {
                    res = 1;
                } else {
                    res = _chain[0](matrix_indices[VDJ_DIV_DEL], matrix_indices[VDJ_CHAIN_SIZE - 1]);
                }
            }
//            res = (*this)(node_i, matrix_i, row_i, column_i);
            return res;
        }


        prob_t backward(node_ind_t node_i,
                        dim_t row_i,
                        dim_t column_i,
                        const vector<matrix_ind_t> &matrix_indices) {
            // find next nodes, compute backward probabilities for them, sum them up and return the sum
            prob_t res = 1;
            node_ind_t next_node_i = node_i + 1;
            if (next_node_i < _chain.size() - 1) {
                for (dim_t next_column_i = 0; next_column_i < nodeColumns(next_node_i); ++next_column_i) {
                    res += (*this)(node_i, matrix_indices[node_i], row_i, column_i) * backward(next_node_i, matrix_indices[next_node_i], column_i, next_column_i);
                }
            } else {
                res = (*this)(node_i, matrix_i, row_i, column_i);
            }
            return res;
        }


        MAAGForwardBackwardAlgorithm() { }

    };

}

#endif //YMIR_MAAGFORWARDBACKWARDALGORITHM_H
