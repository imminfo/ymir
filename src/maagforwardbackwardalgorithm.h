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
            if (_forward_acc) { delete _forward_acc; }
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
        ProbMMC *_forward_acc, *_backward_acc;  /** Temporary MMC for storing forward and backward probabilities correspond. */
        ProbMMC *_fb_acc;  /** Accumulator MMC for storing forward and backward probabilities. */
        prob_t _full_prob;  /** Full generation probability of the input MAAG. */
        size_t _node_i, _mat_i, _row, _column;


        bool init_and_process(const MAAG &maag) {
            _status = false;
            _node_i = 0;
            _mat_i = 0;
            _row = 0;
            _column = 0;
            _full_prob = 0;
            if (maag && maag->_events) {
                _chain = maag._chain;
                _events = maag._events;
                _status = true;
                if (maag.recombination() == VJ_RECOMB) {
                    this->forward_backward_vj(maag);
                } else if (maag.recombination() == VDJ_RECOMB) {
                    this->forward_backward_vdj(maag);
                } else {
                    cerr << "MAAG forward-backward algorithm error: unknown recombination type." << endl;
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
            _forward_acc = new ProbMMC();
            _forward_acc->resize(maag.chainSize());
            // VJ probabilities (for fixed J in future)
            _forward_acc->initNode(VJ_VAR_JOI_GEN_I, 1, maag.nodeRows(VJ_VAR_JOI_GEN_I), 1);
            // V deletions
            _forward_acc->initNode(VJ_VAR_DEL_I, maag.nodeSize(VJ_VAR_DEL_I), maag.nodeRows(VJ_VAR_DEL_I), maag.nodeColumns(VJ_VAR_DEL_I));
            // VJ insertions
            _forward_acc->initNode(VJ_VAR_JOI_INS_I, 1, maag.nodeRows(VJ_VAR_JOI_INS_I), maag.nodeColumns(VJ_VAR_JOI_INS_I));
            // J deletions (for fixed J in future)
            _forward_acc->initNode(VJ_JOI_DEL_I, 1, maag.nodeRows(VJ_JOI_DEL_I), 1);
        }


        // make a matrix chain with forward probabilities for VJ receptors
        void forward_vj(const MAAG &maag, eventind_t j_ind) {
            prob_t temp_prob = 0;
            this->fillZero(_forward_acc);

            // VJ probabilities for the fixed J
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_VAR_JOI_GEN_I); ++row_i) {
                _forward_acc->matrix(VJ_VAR_JOI_GEN_I, 0)(row_i, 0) = maag.matrix(VJ_VAR_JOI_GEN_I, 0)(row_i, j_ind);
            }

            // V deletions
            for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_DEL_I); ++col_i) {
                    _forward_acc->matrix(VJ_VAR_DEL_I, v_ind)(0, col_i) =
                            _forward_acc->matrix(VJ_VAR_JOI_GEN_I, 0)(v_ind, 0) * maag.matrix(VJ_VAR_DEL_I, v_ind)(0, col_i);
                }
            }

            // VD insertions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_VAR_JOI_INS_I); ++row_i) {
                // sum of fi for V del
                for (dim_t v_i = 0; v_i < maag.nodeSize(VJ_VAR_DEL_I); ++v_i) {
                    _forward_acc->matrix(VJ_VAR_JOI_INS_I, 0)(row_i, 0) += _forward_acc->matrix(VJ_VAR_DEL_I, v_i)(0, row_i);
                }
                _forward_acc->matrix(VJ_VAR_JOI_INS_I, 0)(row_i, 0) *= maag.matrix(VJ_VAR_JOI_INS_I, 0)(row_i, 0);

                // for this row in VJ insertions fill with equal values
                // because this rows related to a specific number of deletions for each V gene
                for (dim_t col_i = 1; col_i < maag.nodeColumns(VJ_VAR_JOI_INS_I); ++col_i) {
                    _forward_acc->matrix(VJ_VAR_JOI_INS_I, 0)(row_i, col_i) = _forward_acc->matrix(VJ_VAR_JOI_INS_I, 0)(row_i, 0);
                }
            }

            // J deletions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_JOI_DEL_I); ++row_i) {
                for (dim_t row_vj_i = 0; row_vj_i < maag.nodeSize(VJ_VAR_DEL_I); ++row_vj_i) {
                    _forward_acc->matrix(VJ_JOI_DEL_I, 0)(row_i, 0) += _forward_acc->matrix(VJ_VAR_JOI_INS_I, 0)(row_vj_i, row_i);
                }
                _forward_acc->matrix(VJ_JOI_DEL_I, 0)(row_i, 0) *= maag.matrix(VJ_JOI_DEL_I, j_ind)(row_i, 0);
            }


            // update the full generation probability
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_JOI_DEL_I); ++row_i) {
                _full_prob += _forward_acc->matrix(VJ_JOI_DEL_I, 0)(row_i, 0);
            }
        }


        // make a matrix chain with backward probabilities for VJ receptors
        void backward_vj(const MAAG &maag, eventind_t j_ind) {
            prob_t temp_prob = 0;
            this->fillZero(_backward_acc);

            // J deletions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_JOI_DEL_I); ++row_i) {
                _backward_acc->matrix(VJ_JOI_DEL_I, 0)(row_i, 0) = 1;
            }

            // VJ insertions
            for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_JOI_INS_I); ++col_i) {
                for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_VAR_JOI_INS_I); ++row_i) {
                    _backward_acc->matrix(VJ_VAR_JOI_INS_I)(row_i, col_i) += maag.matrix(VJ_JOI_DEL_I, j_ind)(col_i, 0);
                }
            }

            // V deletions
            for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_DEL_I); ++col_i) {
                    for (dim_t ins_col_i = 0; ins_col_i < maag.nodeColumns(VJ_VAR_JOI_INS_I); ++ins_col_i) {
                        _backward_acc->matrix(VJ_VAR_DEL_I, v_ind)(0, col_i) +=
                                maag.matrix(VJ_VAR_JOI_INS_I, 0)(col_i, ins_col_i) * _backward_acc(VJ_VAR_JOI_INS_I, 0)(col_i, ins_col_i);
                    }
                }
            }

            // V-J genes
            for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_DEL_I); ++col_i) {
                    _backward_acc->matrix(VJ_VAR_JOI_GEN_I, 0)(v_ind, 0) +=
                            _backward_acc->matrix(VJ_VAR_DEL_I, v_ind)(0, col_i) * maag.matrix(VJ_VAR_DEL_I, v_ind)(0, col_i);
                }
            }

        }


        // compute all forward-backward probabilities
        void forward_backward_vj(const MAAG &maag) {
            this->init_vj(maag);

            for (eventind_t j_ind = 0; j_ind <- maag.nJoi(); ++j_ind) {
                // compute forward and backward probabilities for a specific J gene
                this->forward_vj(maag, j_ind);
                this->backward_vj(maag, j_ind);
                // add fi * bi for this J to the accumulator
                for (node_ind_t node_i = 0; node_i < _forward_acc->chainSize(); ++node_i) {
                    for (dim_t row_i = 0; row_i < _forward_acc->matrix(node_i, 0).rows(); ++row_i) {
                        for (dim_t col_i = 0; col_i < _forward_acc->matrix(node_i, 0).cols(); ++col_i) {
                            _fb_acc->matrix(node_i, 0)(row_i, col_i) +=
                                    _forward_acc->matrix(node_i, 0)(row_i, col_i) * _backward_acc->matrix(node_i, 0)(row_i, col_i);
                        }
                    }
                }
            }
        }


        //
        // Forward-backward algorithm for VDJ recombination receptors
        //

        // Initialize accumulator and temporaty MMCs.
        void init_vdj(const MAAG &maag) {
            this->initMMC(maag);

            _forward_acc = new ProbMMC();
        }


        // make a matrix chain with forward probabilities for VDJ receptors
        void forward_vdj(const MAAG &maag) {
            // forward probabilities for (V prob -> V del -> VD ins) are fixed
            // for all pairs of J-D.

            // compute probabilities for every J

            // compute probabilities for every D
        }


        // make a matrix chain with backward probabilities for VDJ receptors
        void backward_vdj(const MAAG &maag) {
            // compute probabilities for every J

            // compute probabilities for every D
        }


        void forward_backward_vdj(const MAAG &maag) {

        }


        MAAGForwardBackwardAlgorithm() { }

    };

}

#endif //YMIR_MAAGFORWARDBACKWARDALGORITHM_H
