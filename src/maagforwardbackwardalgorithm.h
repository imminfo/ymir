//
// Created by Vadim N. on 20/04/2015.
//

#ifndef YMIR_MAAGFORWARDBACKWARDALGORITHM_H
#define YMIR_MAAGFORWARDBACKWARDALGORITHM_H

#include <unordered_map>

#include "maag.h"

namespace ymir {


    class MAAGForwardBackwardAlgorithm {

    protected:

        typedef ProbMMC::dim_t dim_t;
        typedef ProbMMC::matrix_ind_t matrix_ind_t;
        typedef ProbMMC::node_ind_t node_ind_t;

    public:


        MAAGForwardBackwardAlgorithm(const MAAG &maag) {
            init_and_process(maag);
        }


        virtual ~MAAGForwardBackwardAlgorithm() {
            if (_forward_acc) { delete _forward_acc; }
            if (_backward_acc) { delete _backward_acc; }
//            if (_fb_acc) { delete _fb_acc; }
        }


        event_pair_t nextEvent() {
            if (_status) {
                // for each cell of each matrix that is related to
                // an event E return forward(E) * backward(E) / P,
                // where P is the full probability of generation.
            }

            _status = false;
            return event_pair_t(0, 0);
        }


        bool is_empty() {

        }


        bool status() const { return _status; }


        prob_t fullProbability() const { return _full_prob; }


        prob_t bfullProbability() const { return _back_full_prob; }

    protected:

        bool _status;
        ProbMMC *_forward_acc, *_backward_acc;  /** Temporary MMC for storing forward and backward probabilities correspond. */
        prob_t _full_prob;  /** Full generation probability of the input MAAG. */
        prob_t _back_full_prob;  /** Full generation probability of the input MAAG obtained with backward algorithm. Just for testing purposes. */
        vector<event_pair_t> _pairs;


        bool init_and_process(const MAAG &maag) {
            _status = false;
            _pairs.resize(0);
            _full_prob = 0;
            _back_full_prob = 0;
            _forward_acc = nullptr;
            _backward_acc = nullptr;
            _pairs.reserve(maag._chain.size());
            if (maag._events) {
//                _chain = maag._chain;
//                _events = maag._events;
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


        void fillZero(ProbMMC *mmc, uint start_node = 0) {
            for (uint i = start_node; i < mmc->chainSize(); ++i) {
                for (uint j = 0; j < mmc->nodeSize(i); ++j) {
                    mmc->fill(i, j, 0);
                }
            }
        }


        void pushEventPairs(const MAAG &maag, node_ind_t node_i, matrix_ind_t mat_i_fb, matrix_ind_t mat_i_algo) {
            for (dim_t row_i = 0; row_i < maag.nodeRows(node_i); ++row_i) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(node_i); ++col_i) {
                    _pairs.push_back(event_pair_t(
                            maag.event_index(node_i, mat_i_fb, row_i, col_i),
                            (*_forward_acc)(node_i, mat_i_algo, row_i, col_i) * (*_backward_acc)(node_i, mat_i_algo, row_i, col_i)));
                }
            }
        }



        //
        // Forward-backward algorithm for VJ recombination receptors
        //

        // make a matrix chain with forward probabilities for VJ receptors
        void forward_vj(const MAAG &maag, eventind_t j_ind) {
            this->fillZero(_forward_acc);

            // VJ probabilities for the fixed J
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_VAR_JOI_GEN_I); ++row_i) {
                (*_forward_acc)(VJ_VAR_JOI_GEN_I, 0, row_i, 0) = maag(VJ_VAR_JOI_GEN_I, 0, row_i, j_ind);
            }

            // V deletions
            for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_DEL_I); ++col_i) {
                    (*_forward_acc)(VJ_VAR_DEL_I, v_ind, 0, col_i) =
                            (*_forward_acc)(VJ_VAR_JOI_GEN_I, 0, v_ind, 0) * maag(VJ_VAR_DEL_I, v_ind, 0, col_i);
                }
            }

            // VJ insertions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_VAR_JOI_INS_I); ++row_i) {
                prob_t temp_prob = 0;
                // sum of fi for V del
                for (dim_t v_i = 0; v_i < maag.nodeSize(VJ_VAR_DEL_I); ++v_i) {
                    temp_prob += (*_forward_acc)(VJ_VAR_DEL_I, v_i, 0, row_i);
                }

                for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_JOI_INS_I); ++col_i) {
                    (*_forward_acc)(VJ_VAR_JOI_INS_I, 0, row_i, col_i) =
                            temp_prob * maag(VJ_VAR_JOI_INS_I, 0, row_i, col_i);
                }
            }

            // J deletions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_JOI_DEL_I); ++row_i) {
                for (dim_t row_vj_i = 0; row_vj_i < maag.nodeRows(VJ_VAR_JOI_INS_I); ++row_vj_i) {
                    (*_forward_acc)(VJ_JOI_DEL_I, 0, row_i, 0) += (*_forward_acc)(VJ_VAR_JOI_INS_I, 0, row_vj_i, row_i);
                }
                (*_forward_acc)(VJ_JOI_DEL_I, 0, row_i, 0) *= maag(VJ_JOI_DEL_I, j_ind, row_i, 0);
            }


            // update the full generation probability
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_JOI_DEL_I); ++row_i) {
                _full_prob += (*_forward_acc)(VJ_JOI_DEL_I, 0, row_i, 0);
            }
        }


        // make a matrix chain with backward probabilities for VJ receptors
        void backward_vj(const MAAG &maag, eventind_t j_ind) {
            this->fillZero(_backward_acc);

            // J deletions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_JOI_DEL_I); ++row_i) {
                (*_backward_acc)(VJ_JOI_DEL_I, 0, row_i, 0) = 1;
            }

            // VJ insertions
            for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_JOI_INS_I); ++col_i) {
                for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_VAR_JOI_INS_I); ++row_i) {
                    (*_backward_acc)(VJ_VAR_JOI_INS_I, 0, row_i, col_i) += maag(VJ_JOI_DEL_I, j_ind, col_i, 0);
                }
            }

            // V deletions
            for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_DEL_I); ++col_i) {
                    for (dim_t ins_col_i = 0; ins_col_i < maag.nodeColumns(VJ_VAR_JOI_INS_I); ++ins_col_i) {
                        (*_backward_acc)(VJ_VAR_DEL_I, v_ind, 0, col_i) +=
                                (*_backward_acc)(VJ_VAR_JOI_INS_I, 0, col_i, ins_col_i) * maag(VJ_VAR_JOI_INS_I, 0, col_i, ins_col_i);
                    }
                }
            }

            // V-J genes
            for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_DEL_I); ++col_i) {
                    (*_backward_acc)(VJ_VAR_JOI_GEN_I, 0, v_ind, 0) +=
                            (*_backward_acc)(VJ_VAR_DEL_I, v_ind, 0, col_i) * maag(VJ_VAR_DEL_I, v_ind, 0, col_i);
                }
            }

            for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                _back_full_prob +=
                        (*_backward_acc)(VJ_VAR_JOI_GEN_I, 0, v_ind, 0) * maag(VJ_VAR_JOI_GEN_I, 0, v_ind, j_ind);
            }
        }


        // compute all forward-backward probabilities
        void forward_backward_vj(const MAAG &maag) {
            _forward_acc = new ProbMMC();
            _forward_acc->resize(maag.chainSize());
            // VJ probabilities (for fixed J in future)
            _forward_acc->initNode(VJ_VAR_JOI_GEN_I, 1, maag.nodeRows(VJ_VAR_JOI_GEN_I), 1);
            // V deletions
            _forward_acc->initNode(VJ_VAR_DEL_I, maag.nodeSize(VJ_VAR_DEL_I), 1, maag.nodeColumns(VJ_VAR_DEL_I));
            // VJ insertions
            _forward_acc->initNode(VJ_VAR_JOI_INS_I, 1, maag.nodeRows(VJ_VAR_JOI_INS_I), maag.nodeColumns(VJ_VAR_JOI_INS_I));
            // J deletions (for fixed J in future)
            _forward_acc->initNode(VJ_JOI_DEL_I, 1, maag.nodeRows(VJ_JOI_DEL_I), 1);

            _backward_acc = new ProbMMC();
            _backward_acc->resize(maag.chainSize());
            // VJ probabilities (for fixed J in future)
            _backward_acc->initNode(VJ_VAR_JOI_GEN_I, 1, maag.nodeRows(VJ_VAR_JOI_GEN_I), 1);
            // V deletions
            _backward_acc->initNode(VJ_VAR_DEL_I, maag.nodeSize(VJ_VAR_DEL_I), 1, maag.nodeColumns(VJ_VAR_DEL_I));
            // VJ insertions
            _backward_acc->initNode(VJ_VAR_JOI_INS_I, 1, maag.nodeRows(VJ_VAR_JOI_INS_I), maag.nodeColumns(VJ_VAR_JOI_INS_I));
            // J deletions (for fixed J in future)
            _backward_acc->initNode(VJ_JOI_DEL_I, 1, maag.nodeRows(VJ_JOI_DEL_I), 1);

            // Compute fi * bi / Pgen for each event.
            for (eventind_t j_ind = 0; j_ind < maag.nJoi(); ++j_ind) {
                // compute forward and backward probabilities for a specific J gene
                this->forward_vj(maag, j_ind);
                this->backward_vj(maag, j_ind);

                // add fi * bi for this J to the accumulator
                for (node_ind_t node_i = 0; node_i < _forward_acc->chainSize(); ++node_i) {
                    for (dim_t row_i = 0; row_i < _forward_acc->nodeRows(node_i); ++row_i) {
                        for (dim_t col_i = 0; col_i < _forward_acc->nodeColumns(node_i); ++col_i) {
//                            (*_fb_acc)(node_i, 0, row_i, col_i) +=
//                                    (*_forward_acc)(node_i, 0, row_i, col_i) * (*_backward_acc)(node_i, 0, row_i, col_i);
                        }
                    }
                }
            }
        }


        //
        // Forward-backward algorithm for VDJ recombination receptors
        //

        // make a matrix chain with forward probabilities for VDJ receptors
        void forward_vdj(const MAAG &maag, eventind_t d_ind, eventind_t j_ind, bool recompute_d_gen_fi) {
            if (recompute_d_gen_fi) {
                // forward probabilities for (V prob -> V del -> VD ins) are fixed
                // for all pairs of J-D.
                // We have already stored in _forward_acc fi for this D, so we don't
                // need to recompute entire _forward_acc, we just need to recompute
                // J deletions and J genes fi.
                this->fillZero(_forward_acc, VDJ_DIV_DEL_I);

                // D deletions
                for (dim_t row_i = 0; row_i < maag.nodeRows(VDJ_DIV_DEL_I); ++row_i) {
                    for (dim_t col_i = 0; col_i < maag.nodeColumns(VDJ_DIV_DEL_I); ++col_i) {
                        for (dim_t ins_row_i = 0; ins_row_i < maag.nodeRows(VDJ_VAR_DIV_INS_I); ++ins_row_i) {
                            (*_forward_acc)(VDJ_DIV_DEL_I, 0, row_i, col_i) +=
                                    (*_forward_acc)(VDJ_VAR_DIV_INS_I, 0, ins_row_i, row_i) * maag(VDJ_DIV_DEL_I, d_ind, row_i, col_i);
                        }
                    }
                }

                // DJ insertions
                for (dim_t row_i = 0; row_i < maag.nodeRows(VDJ_DIV_JOI_INS_I); ++row_i) {
                    for (dim_t col_i = 0; col_i < maag.nodeColumns(VDJ_DIV_JOI_INS_I); ++col_i) {
                        for (dim_t dgen_row_i = 0; dgen_row_i < maag.nodeRows(VDJ_DIV_DEL_I); ++dgen_row_i) {
                            (*_forward_acc)(VDJ_DIV_JOI_INS_I, 0, row_i, col_i) +=
                                    (*_forward_acc)(VDJ_DIV_DEL_I, 0, dgen_row_i, row_i) * maag(VDJ_DIV_JOI_INS_I, 0, row_i, col_i);
                        }
                    }
                }

            } else {
                this->fillZero(_forward_acc, VDJ_JOI_DEL_I);
            }

            // J deletions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VDJ_JOI_DEL_I); ++row_i) {
                for (dim_t ins_row_i = 0; ins_row_i < maag.nodeRows(VDJ_DIV_JOI_INS_I); ++ins_row_i) {
                    (*_forward_acc)(VDJ_JOI_DEL_I, 0, row_i, 0) +=
                            (*_forward_acc)(VDJ_DIV_JOI_INS_I, 0, ins_row_i, row_i) * maag(VDJ_JOI_DEL_I, j_ind, row_i, 0);
                }
            }

            // J-D genes
            for (dim_t row_i = 0; row_i < maag.nodeRows(VDJ_JOI_DEL_I); ++row_i) {
                (*_forward_acc)(VDJ_JOI_DIV_GEN_I, 0, 0, 0) +=
                        (*_forward_acc)(VDJ_JOI_DEL_I, 0, row_i, 0) * maag(VDJ_JOI_DIV_GEN_I, 0, j_ind, d_ind);
            }
//            (*_forward_acc)(VDJ_JOI_DIV_GEN_I)(0, 0) =
//                    (*_forward_acc)(VDJ_JOI_DEL_I).sum() * maag.matrix(VDJ_JOI_DIV_GEN_I)(j_ind, d_ind);

            // update the full generation probability
            _full_prob += (*_forward_acc)(VDJ_JOI_DIV_GEN_I, 0, 0, 0);
        }


        // make a matrix chain with backward probabilities for VDJ receptors
        void backward_vdj(const MAAG &maag, eventind_t d_ind, eventind_t j_ind) {
            this->fillZero(_backward_acc);

            // J-D pairs
            (*_backward_acc)(VDJ_JOI_DIV_GEN_I, 0, 0, 0) = 1;

            // J deletions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VDJ_JOI_DEL_I); ++row_i) {
                (*_backward_acc)(VDJ_JOI_DEL_I, 0, row_i, 0) += maag(VDJ_JOI_DIV_GEN_I, 0, j_ind, d_ind);
            }

            // DJ insertions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VDJ_DIV_JOI_INS_I); ++row_i) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VDJ_DIV_JOI_INS_I); ++col_i) {
                    (*_backward_acc)(VDJ_DIV_JOI_INS_I, 0, row_i, col_i) +=
                            (*_backward_acc)(VDJ_JOI_DEL_I, 0, col_i, 0) * maag(VDJ_JOI_DEL_I, j_ind, col_i, 0);
                }
            }

            // D5'-3' deletions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VDJ_DIV_DEL_I); ++row_i) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VDJ_DIV_DEL_I); ++col_i) {
                    for (dim_t ins_col_i = 0; ins_col_i < maag.nodeColumns(VDJ_DIV_JOI_INS_I); ++ins_col_i) {
                        (*_backward_acc)(VDJ_DIV_DEL_I, 0, row_i, col_i) +=
                                (*_backward_acc)(VDJ_DIV_JOI_INS_I, 0, col_i, ins_col_i) * maag(VDJ_DIV_JOI_INS_I, 0, col_i, ins_col_i);
                    }
                }
            }

            // VD insertions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VDJ_VAR_DIV_INS_I); ++row_i) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VDJ_VAR_DIV_INS_I); ++col_i) {
                    for (dim_t d_col_i = 0; d_col_i < maag.nodeColumns(VDJ_DIV_DEL_I); ++d_col_i) {
                        (*_backward_acc)(VDJ_VAR_DIV_INS_I, 0, row_i, col_i) +=
                                (*_backward_acc)(VDJ_DIV_DEL_I, 0, col_i, d_col_i) * maag(VDJ_DIV_DEL_I, d_ind, col_i, d_col_i);
                    }
                }
            }

            // V deletions and V genes
            for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VDJ_VAR_DEL_I); ++col_i) {
                    for (dim_t ins_col_i = 0; ins_col_i < maag.nodeColumns(VDJ_VAR_DIV_INS_I); ++ins_col_i) {
                        (*_backward_acc)(VDJ_VAR_DEL_I, v_ind, 0, col_i) +=
                                (*_backward_acc)(VDJ_VAR_DIV_INS_I, 0, col_i, ins_col_i) * maag(VDJ_VAR_DIV_INS_I, 0, col_i, ins_col_i);
                    }
                    (*_backward_acc)(VDJ_VAR_GEN_I, v_ind, 0, 0) +=
                            (*_backward_acc)(VDJ_VAR_DEL_I, v_ind, 0, col_i) * maag(VDJ_VAR_DEL_I, v_ind, 0, col_i);
                }

                // update the full (back) generation probability
                _back_full_prob += (*_backward_acc)(VDJ_VAR_GEN_I, v_ind, 0, 0) * maag(VDJ_VAR_GEN_I, v_ind, 0, 0);
            }
        }


        void forward_backward_vdj(const MAAG &maag) {
            _forward_acc = new ProbMMC();
            _forward_acc->resize(maag.chainSize());
            // V genes - fi which is constant for all J-D pairs
            _forward_acc->initNode(VDJ_VAR_GEN_I, maag.nodeSize(VDJ_VAR_GEN_I), 1, 1);
            // V deletions - fi which is constant for all J-D pairs
            _forward_acc->initNode(VDJ_VAR_DEL_I, maag.nodeSize(VDJ_VAR_DEL_I), 1, maag.nodeColumns(VDJ_VAR_DEL_I));
            // VD insertions
            _forward_acc->initNode(VDJ_VAR_DIV_INS_I, 1, maag.nodeRows(VDJ_VAR_DIV_INS_I), maag.nodeColumns(VDJ_VAR_DIV_INS_I));
            // D5' - 3' deletions
            _forward_acc->initNode(VDJ_DIV_DEL_I, 1, maag.nodeRows(VDJ_DIV_DEL_I), maag.nodeColumns(VDJ_DIV_DEL_I));
            // DJ insertions
            _forward_acc->initNode(VDJ_DIV_JOI_INS_I, 1, maag.nodeRows(VDJ_DIV_JOI_INS_I), maag.nodeColumns(VDJ_DIV_JOI_INS_I));
            // J deletions
            _forward_acc->initNode(VDJ_JOI_DEL_I, 1, maag.nodeRows(VDJ_JOI_DEL_I), maag.nodeColumns(VDJ_JOI_DEL_I));
            // J-D pairs
//            _forward_acc->initNode(VDJ_JOI_DIV_GEN_I, 1, maag.nodeRows(VDJ_JOI_DIV_GEN_I), 1);
            _forward_acc->initNode(VDJ_JOI_DIV_GEN_I, 1, 1, 1);

            _backward_acc = new ProbMMC();
            _backward_acc->resize(maag.chainSize());
            // V genes
            _backward_acc->initNode(VDJ_VAR_GEN_I, maag.nodeSize(VDJ_VAR_GEN_I), 1, 1);
            // V deletions
            _backward_acc->initNode(VDJ_VAR_DEL_I, maag.nodeSize(VDJ_VAR_DEL_I), 1, maag.nodeColumns(VDJ_VAR_DEL_I));
            // VD insertions
            _backward_acc->initNode(VDJ_VAR_DIV_INS_I, 1, maag.nodeRows(VDJ_VAR_DIV_INS_I), maag.nodeColumns(VDJ_VAR_DIV_INS_I));
            // D5' - 3' deletions
            _backward_acc->initNode(VDJ_DIV_DEL_I, 1, maag.nodeRows(VDJ_DIV_DEL_I), maag.nodeColumns(VDJ_DIV_DEL_I));
            // DJ insertions
            _backward_acc->initNode(VDJ_DIV_JOI_INS_I, 1, maag.nodeRows(VDJ_DIV_JOI_INS_I), maag.nodeColumns(VDJ_DIV_JOI_INS_I));
            // J deletions
            _backward_acc->initNode(VDJ_JOI_DEL_I, 1, maag.nodeRows(VDJ_JOI_DEL_I), maag.nodeColumns(VDJ_JOI_DEL_I));
            // J-D pairs
//            _backward_acc->initNode(VDJ_JOI_DIV_GEN_I, 1, maag.nodeRows(VDJ_JOI_DIV_GEN_I), 1);
            _backward_acc->initNode(VDJ_JOI_DIV_GEN_I, 1, 1, 1);

            // Because fi for V genes, V deletions and VD insertions are constant for all
            // pairs of J-D, we compute them here.
            this->fillZero(_forward_acc);
            // V genes and deletions
            for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                // gene probability
                (*_forward_acc)(VDJ_VAR_GEN_I, v_ind, 0, 0) = maag(VDJ_VAR_GEN_I, v_ind, 0, 0);

                // deletions probabilities
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VDJ_VAR_DEL_I); ++col_i) {
                    (*_forward_acc)(VDJ_VAR_DEL_I, v_ind, 0, col_i) =
                            (*_forward_acc)(VDJ_VAR_GEN_I, v_ind, 0, 0) * maag(VDJ_VAR_DEL_I, v_ind, 0, col_i);
                }
            }
            // VD insertions
            for (dim_t row_i = 0; row_i < maag.nodeRows(VDJ_VAR_DIV_INS_I); ++row_i) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VDJ_VAR_DIV_INS_I); ++col_i) {
                    for (eventind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                        (*_forward_acc)(VDJ_VAR_DIV_INS_I, 0, row_i, col_i) +=
                                (*_forward_acc)(VDJ_VAR_DEL_I, v_ind, 0, row_i) * maag(VDJ_VAR_DIV_INS_I, 0, row_i, col_i);
                    }
                }
            }

            // Compute fi * bi / Pgen for each event.
            bool recompute_d_gen_fi = true;
            for (eventind_t d_ind = 0; d_ind < maag.nDiv(); ++d_ind) {
                recompute_d_gen_fi = true;
                for (eventind_t j_ind = 0; j_ind < maag.nJoi(); ++j_ind) {
                    // compute forward and backward probabilities for a specific J gene
                    this->forward_vdj(maag, d_ind, j_ind, recompute_d_gen_fi);
                    this->backward_vdj(maag, d_ind, j_ind);
                    recompute_d_gen_fi = false;

                    // add fi * bi for this J to the accumulator
                    for (node_ind_t node_i = 0; node_i < _forward_acc->chainSize(); ++node_i) {
                        for (dim_t row_i = 0; row_i < _forward_acc->nodeRows(node_i); ++row_i) {
                            for (dim_t col_i = 0; col_i < _forward_acc->nodeColumns(node_i); ++col_i) {
//                                (*_fb_acc)(node_i, 0, row_i, col_i) +=
//                                        (*_forward_acc)(node_i, 0, row_i, col_i) * (*_backward_acc)(node_i, 0, row_i, col_i);
                            }
                        }
                    }
                }
            }
        }


        MAAGForwardBackwardAlgorithm() { }

    };

}

#endif //YMIR_MAAGFORWARDBACKWARDALGORITHM_H
