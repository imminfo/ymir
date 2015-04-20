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
                _events = maag._events;
                _status = true;
            }
            return _status;
        }


        event_pair_t next() {
            if (_status && _node_i < _events->chainSize()) {
                // node indices
                if (_mat_i >= _events->nodeSize()) {
                    ++_node_i;
                    _mat_i = 0;
                    _row = 0;
                    _column = 0;
                    return next();
                } else {

                }

                // matrix indices
                if (_row >= _events->matrix(_node_i, _mat_i).rows()) {
                    ++_mat_i;
                    _row = 0;
                    _column = 0;
                    return next();
                } else {

                }

                // column indices
                if (_column >= _events->matrix(_node_i, _mat_i).cols()) {
                    ++_row;
                    _column = 0;
                    return next();
                } else {

                    ++_column;
                }
            }

            _status = false;
            retrun event_pair_t(0, 0);
        }


        bool status() const { return _status; }

    protected:

        bool _status;
        size_t _node_i, _mat_i, _row, _column;


        MAAGForwardBackwardAlgorithm() { }

    };

}

#endif //YMIR_MAAGFORWARDBACKWARDALGORITHM_H
