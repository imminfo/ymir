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
            init_and_process(maag);
        }


        virtual ~MAAGForwardBackwardAlgorithm() {
            if (_temp_mmc) { delete _temp_mmc; }
            if (_forward_acc) { delete _forward_acc; }
            if (_backward_acc) { delete _backward_acc; }
        }


        event_pair_t next() {
            if (_status) {
                // for each cell of each matrix that is related to
                // an event E return forward(E) * backward(E) / P,
                // where P is the full probability of generation.
            }

            _status = false;
            return event_pair_t(0, 0);
        }


        bool status() const { return _status; }

    protected:

        bool _status;
        bool _processed;
        RECOMBINATION _recomb;
        size_t _node_i, _mat_i, _row, _column;
        ProbMMC *_temp_mmc;  /** Temporary MMC for storing temporary results for forward and backward probabilities */
        ProbMMC *_forward_acc, *_backward_acc;  /** Accumulator MMC for forward and backward probabilities */


        bool init_and_process(const MAAG &maag) {
            _status = false;
            _processed = false;
            _node_i = 0;
            _mat_i = 0;
            _row = 0;
            _column = 0;
            if (maag && maag->_events) {
                _chain = maag._chain;
                _events = maag._events;
                _status = true;
                _recomb = maag.recombination();
                if (_recomb == VJ_RECOMB) {
                    this->init_vj(maag);
                    // Should I make it lazy?
                    this->forward_vj(maag);
                    this->backward_vj(maag);
                } else if (_recomb == VDJ_RECOMB) {
                    this->init_vdj(maag);
                    this->forward_vdj(maag);
                    this->backward_vdj(maag);
                } else {
                    cerr << "Forward-backward algorithm error: unknown recombination type." << endl;
                    _status = false;
                }
            }
            return _status;
        }


        void copyMaagShape(const MAAG &source, ProbMMC *target) {
            target = new ProbMMC();
            target->resize(source.chainSize());
            for (uint i = 0; i < source.chainSize(); ++i) {
                target->initNode(i, source.nodeSize(i), source.nodeRows(i), source.nodeColumns(i));
            }
        }


        void initMMC(const MAAG &maag) {
            copyMaagShape(maag, _forward_acc);
            copyMaagShape(maag, _backward_acc);
        }


        void fillZero(ProbMMC *mmc) {
            for (uint i = 0; i < mmc->chainSize(); ++i) {
                for (uint j = 0; j < mmc->nodeSize(i); ++j) {
                    mmc->matrix(i, j).fill(0);
                }
            }
        }


        //
        // Forward-backward algorithm for VJ recombination receptors
        //

        // Initialize accumulator and temporaty MMCs.
        void init_vj(const MAAG &maag) {
            this->initMMC(maag);

            // Initialise the temporary MMC for computing forward-backward probabilities.
            _temp_mmc = new ProbMMC();
            _temp_mmc->resize(maag.chainSize());
            // VJ probabilities (for fixed J in future)
            _temp_mmc->initNode(0, 1, maag.nodeRows(0), 1);
            // V deletions
            _temp_mmc->initNode(1, maag.nodeSize(1), maag.nodeRows(1), maag.nodeColumns(1));
            // VJ insertions
            _temp_mmc->initNode(2, 1, maag.nodeRows(2), maag.nodeColumns(2));
            // J deletions (for fixed J in future)
            _temp_mmc->initNode(3, 1, maag.nodeRows(3), 1);
        }


        // make a matrix chain with forward probabilities for VJ receptors
        prob_t forward_vj(const MAAG &maag) {
            // and compute probabilities for every J
            prob_t temp_prob = 0;
            for (eventind_t j_ind = 0; j_ind < maag.nJoi(); ++j_ind) {
                this->fillZero(_temp_mmc);

                // VJ probabilities for the fixed J
                for (dim_t row_i = 0; row_i < _temp_mmc->nodeRows(0); ++row_i) {
                    _temp_mmc->matrix(0, 0)(row_i, 0) = maag.matrix(0, 0)(row_i, j_ind);
                }

                // V deletions
                for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                    for (dim_t col_i = 0; col_i < _temp_mmc->nodeColumns(1); ++col_i) {
                        _temp_mmc->matrix(1, v_ind)(0, col_i) = maag.matrix(0, 0)(v_ind, j_ind) * maag.matrix(0, v_ind)(0, col_i);
                    }
                }

                // VD insertions
                for (dim_t row_i = 0; row_i < _temp_mmc->nodeRows(2); ++row_i) {
                    for (dim_t col_i = 0; col_i < _temp_mmc->nodeColumns(2); ++col_i) {
                        // sum of fi for V del
                        for (dim_t prev_v_col_i = 0;)
                        _temp_mmc->matrix(2, 0)(row_i, col_i) = ;
                    }
                }

                // J deletions

            }
        }


        // make a matrix chain with backward probabilities for VJ receptors
        prob_t backward_vj(const MAAG &maag) {
            // compute probabilities for every J
        }


        //
        // Forward-backward algorithm for VDJ recombination receptors
        //

        // Initialize accumulator and temporaty MMCs.
        void init_vdj(const MAAG &maag) {
            this->initMMC(maag);

            _temp_mmc = new ProbMMC();
        }


        // make a matrix chain with forward probabilities for VDJ receptors
        prob_t forward_vdj(const MAAG &maag) {
            // forward probabilities for (V prob -> V del -> VD ins) are fixed
            // for all pairs of J-D.

            // compute probabilities for every J

            // compute probabilities for every D
        }


        // make a matrix chain with backward probabilities for VDJ receptors
        prob_t backward_vdj(const MAAG &maag) {
            // compute probabilities for every J

            // compute probabilities for every D
        }


        MAAGForwardBackwardAlgorithm() { }

    };

}

#endif //YMIR_MAAGFORWARDBACKWARDALGORITHM_H
