//
// Created by Vadim N. on 20/04/2015.
//

#ifndef YMIR_MAAGFORWARDBACKWARDALGORITHM_H
#define YMIR_MAAGFORWARDBACKWARDALGORITHM_H

#include <unordered_map>

#include "maag.h"

namespace ymir {

    class MAAGForwardBackwardAlgorithm;
    class ForwardBackwardAlgorithm;
    class VJForwardBackward;
    class VDJForwardBackward;


    class MAAGForwardBackwardAlgorithm {

    protected:

        typedef ProbMMC::dim_t dim_t;
        typedef ProbMMC::matrix_ind_t matrix_ind_t;
        typedef ProbMMC::node_ind_t node_ind_t;

    public:


        MAAGForwardBackwardAlgorithm()
            : _pairs_i(0),
              _status(false),
              _vectorised(false),
              _err_prob(0)
        {
        }


        MAAGForwardBackwardAlgorithm(const MAAG &maag, ErrorMode error_mode = NO_ERRORS) {
            process(maag, error_mode);
        }


        virtual ~MAAGForwardBackwardAlgorithm()
        {
        }


        /**
         *
         */
        bool process(const MAAG &maag, ErrorMode error_mode = NO_ERRORS) {
            _pairs_i = 0;
            _status = false;
            _vectorised = false;
            _pairs.clear();
            _pairs.reserve(maag._chain.size() + 10);
            _full_prob = 0;
            _back_full_prob = 0;
            _err_prob = 0;
            _err_mode = error_mode;

            if (_err_mode == COMPUTE_ERRORS && !maag.has_errors()) {
                cerr << "MAAG forward-backward algorithm error: no error matrix has been found in the input MAAG." << endl;
                _status = false;
                return false;
            } else {
                if (maag.has_events()) {
                    _status = true;
                    if (maag.recombination() == VJ_RECOMB) {
                        fill(_nuc_arr1, _nuc_arr1 + 4, 0);
                        this->forward_backward_vj(maag);
                    } else if (maag.recombination() == VDJ_RECOMB) {
                        fill(_nuc_arr1, _nuc_arr1 + 16, 0);
                        fill(_nuc_arr2, _nuc_arr2 + 16, 0);
                        this->forward_backward_vdj(maag);
                    } else {
                        cerr << "MAAG forward-backward algorithm error: unknown recombination type." << endl;
                        _status = false;
                    }
                } else {
                    cerr << "MAAG forward-backward algorithm error: no event matrix has been found in the input MAAG." << endl;
                    _status = false;
                }
            }

            return _status;
        }


        event_pair_t nextEvent() {
            if (_status) {
                event_pair_t res = _pairs[_pairs_i];
//                cout << res.second << " -> ";
//                res.second = res.second / _full_prob;
//                cout << res.second << endl;
                ++_pairs_i;
                return res;
            }

            _status = false;
            return event_pair_t(0, 0);
        }


        bool is_empty() {
            if (_status) {
                if (_pairs_i != _pairs.size()) {
                    return false;
                }
            }
            return true;
        }


        bool status() const { return _status; }


        prob_t fullProbability() const { return _full_prob; }


        prob_t bfullProbability() const { return _back_full_prob; }


        const vector<event_pair_t>& event_pairs() const { return _pairs; }


        prob_t* VJ_nuc_probs() { return _nuc_arr1; }
        prob_t* VD_nuc_probs() { return _nuc_arr1; }
        prob_t* DJ_nuc_probs() { return _nuc_arr2; }

        prob_t err_prob() { return _err_prob; }

    protected:

        bool _status;
        bool _vectorised;
        pProbMMC _forward_acc, _backward_acc;  /** Temporary MMC for storing forward and backward probabilities correspond. */
        prob_t _full_prob;  /** Full generation probability of the input MAAG. */
        prob_t _back_full_prob;  /** Full generation probability of the input MAAG obtained with backward algorithm. Just for testing purposes. */

        vector<event_pair_t> _pairs;
        size_t _pairs_i;
        unordered_map<event_ind_t, prob_t> _pair_map;
        prob_t _nuc_arr1[16], _nuc_arr2[16];
        prob_t _err_prob;
        ErrorMode _err_mode;


        /**
         *
         */
        void fillZero(ProbMMC *mmc, uint start_node = 0) {
            for (uint i = start_node; i < mmc->chainSize(); ++i) {
                for (uint j = 0; j < mmc->nodeSize(i); ++j) {
                    mmc->fill(i, j, 0);
                }
            }
        }


        /**
         *
         */
        void inferInsertionNucleotides(const MAAG &maag, node_ind_t ins_node,
                                       seq_len_t left_start_pos, seq_len_t left_end_pos,
                                       seq_len_t right_start_pos, seq_len_t right_end_pos,
                                       prob_t *nuc_arr, bool reversed = false) {
            node_ind_t forw_node = ins_node - 1, back_node = ins_node;
            prob_t scenario_prob = 0;
            prob_t temp_arr[16];
            uint n = 0;
            seq_len_t left_pos, right_pos, start_shift;

//#ifdef YDEBUG
//            if (left_end_pos - left_start_pos + 1 != maag.nodeRows(forw_node)) { throw(std::runtime_error("Wrong position boundaries (forward node)!")); }
//            if (right_end_pos - right_start_pos + 1 != maag.nodeColumns(back_node)) { throw(std::runtime_error("Wrong position boundaries (backward node)!")); }
//#endif

            if (maag.is_vj()) {
                for (dim_t row_i = 0; row_i < maag.nodeRows(ins_node); ++row_i) {
                    left_pos = left_start_pos + row_i;
                    for (dim_t col_i = 0; col_i < maag.nodeColumns(ins_node); ++col_i) {
                        right_pos = right_start_pos + col_i;

                        if (maag.position(right_pos) - maag.position(left_pos) - 1 > 0 && maag.event_index(ins_node, 0, row_i, col_i)) {
                            fill(temp_arr, temp_arr + 4, 0);
                            n = 0;

                            if (_err_mode == NO_ERRORS) {
                                for (seq_len_t pos = maag.position(left_pos) + 1; pos < maag.position(right_pos); ++pos) {
                                    temp_arr[nuc_hash(maag.sequence()[pos - 1])] += 1;
                                    ++n;
                                }
                            } else {
                                for (seq_len_t pos = maag.position(left_pos) + 1; pos < maag.position(right_pos); ++pos) {
                                    temp_arr[0] += _err_prob / 3;
                                    temp_arr[1] += _err_prob / 3;
                                    temp_arr[2] += _err_prob / 3;
                                    temp_arr[3] += _err_prob / 3;
                                    temp_arr[nuc_hash(maag.sequence()[pos - 1])] += (1 - _err_prob / 3);
                                    ++n;
                                }
                            }

                            scenario_prob = (*_forward_acc)(ins_node, 0, row_i, col_i)
                                            * (*_backward_acc)(back_node, 0, row_i, col_i);

                            for (auto i = 0; i < 4; ++i) {
                                nuc_arr[i] += (temp_arr[i] * scenario_prob) / n;
                            }

//                            _err_prob += scenario_prob / (maag.position(right_pos) - maag.position(left_pos) - 1);

//                            _nuc_arr1[0] /= _cur_prob;
//                            _nuc_arr1[1] /= _cur_prob;
//                            _nuc_arr1[2] /= _cur_prob;
//                            _nuc_arr1[3] /= _cur_prob;
                        }
                    }
                }
            } else {
                for (dim_t row_i = 0; row_i < maag.nodeRows(ins_node); ++row_i) {
                    left_pos = left_start_pos + row_i;
                    for (dim_t col_i = 0; col_i < maag.nodeColumns(ins_node); ++col_i) {
                        right_pos = right_start_pos + col_i;

                        if (maag.position(right_pos) - maag.position(left_pos) - 1 > 0) {
                            fill(temp_arr, temp_arr + 16, 0);
                            n = 0;
                            start_shift = 0;

                            if (!reversed) {
                                if (maag.position(left_pos) == 0) {
                                    start_shift = 1;
                                    temp_arr[4 * nuc_hash('A') + nuc_hash(maag.sequence()[0])] = .25;
                                    temp_arr[4 * nuc_hash('C') + nuc_hash(maag.sequence()[0])] = .25;
                                    temp_arr[4 * nuc_hash('G') + nuc_hash(maag.sequence()[0])] = .25;
                                    temp_arr[4 * nuc_hash('T') + nuc_hash(maag.sequence()[0])] = .25;
                                    ++n;
                                }

                                if (_err_mode == NO_ERRORS) {
                                    for (seq_len_t pos = maag.position(left_pos) + start_shift + 1; pos < maag.position(right_pos); ++pos) {
                                        temp_arr[4 * nuc_hash(maag.sequence()[pos - 2]) + nuc_hash(maag.sequence()[pos - 1])] += 1;
                                        ++n;
                                    }
                                } else {
                                    // TODO: ???
                                }
                            } else {
                                if (maag.position(right_pos) == maag.sequence().size() + 1) {
                                    start_shift = 1;
                                    temp_arr[4 * nuc_hash('A') + nuc_hash(maag.sequence()[maag.sequence().size() - 1])] = .25;
                                    temp_arr[4 * nuc_hash('C') + nuc_hash(maag.sequence()[maag.sequence().size() - 1])] = .25;
                                    temp_arr[4 * nuc_hash('G') + nuc_hash(maag.sequence()[maag.sequence().size() - 1])] = .25;
                                    temp_arr[4 * nuc_hash('T') + nuc_hash(maag.sequence()[maag.sequence().size() - 1])] = .25;
                                    ++n;
                                }

                                for (seq_len_t pos = maag.position(right_pos) - start_shift; pos > maag.position(left_pos) + 1; --pos) {
                                    temp_arr[4 * nuc_hash(maag.sequence()[pos - 1]) + nuc_hash(maag.sequence()[pos - 2])] += 1;
                                    ++n;
                                }
                            }

                            scenario_prob = (*_forward_acc)(ins_node, 0, row_i, col_i)
                                            * (*_backward_acc)(back_node, 0, row_i, col_i);

                            for (auto i = 0; i < 16; ++i) {
                                nuc_arr[i] += (temp_arr[i] * scenario_prob) / n;
                            }
                        }
                    }
                }
            }
        }


        /**
         * \brief Access to a hash map which maps event probabilities to event indices.
         */
        ///@{
        void pushEventValue(event_ind_t event_index, prob_t prob_value) {
            if (prob_value && event_index) {

                auto elem = _pair_map.find(event_index);

                if (elem == _pair_map.end()) {
                    _pair_map[event_index] = 0;
                }
                _pair_map[event_index] += prob_value;
            }
        }

        void pushEventPair(const MAAG &maag, node_ind_t node_i, matrix_ind_t maag_mat_i, dim_t maag_row_i, dim_t maag_col_i,
                                                                matrix_ind_t fb_mat_i, dim_t fb_row_i, dim_t fb_col_i) {
            this->pushEventValue(maag.event_index(node_i, maag_mat_i, maag_row_i, maag_col_i),
                                 (*_forward_acc)(node_i, fb_mat_i, fb_row_i, fb_col_i) * (*_backward_acc)(node_i, fb_mat_i, fb_row_i, fb_col_i));
        }

        void pushEventPairs(const MAAG &maag, node_ind_t node_i, matrix_ind_t maag_mat_i, matrix_ind_t fb_mat_i) {
            for (dim_t row_i = 0; row_i < maag.nodeRows(node_i); ++row_i) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(node_i); ++col_i) {
                    this->pushEventPair(maag, node_i, maag_mat_i, row_i, col_i, fb_mat_i, row_i, col_i);
                }
            }
        }

        void pushEventPairsWithErrors(const MAAG &maag, node_ind_t node_i, matrix_ind_t maag_mat_i, matrix_ind_t fb_mat_i, node_ind_t err_node_i) {
            this->pushEventPairs(maag, node_i, maag_mat_i, fb_mat_i);

            if (_err_mode == COMPUTE_ERRORS) {
                for (dim_t row_i = 0; row_i < maag.nodeRows(node_i); ++row_i) {
                    for (dim_t col_i = 0; col_i < maag.nodeColumns(node_i); ++col_i) {
                        if (maag.errors(err_node_i, fb_mat_i, row_i, col_i)) {
                            _err_prob += (*_forward_acc)(node_i, fb_mat_i, row_i, col_i)
                                         * (*_backward_acc)(node_i, fb_mat_i, row_i, col_i)
                                         / maag.errors(err_node_i, fb_mat_i, row_i, col_i);
                        }
                    }
                }
            }
        }
        ///@}


        /**
         *
         */
        void vectorise_pair_map(const MAAG &maag) {
            if (maag.is_vj()) {
                _nuc_arr1[0] /= _full_prob;
                _nuc_arr1[1] /= _full_prob;
                _nuc_arr1[2] /= _full_prob;
                _nuc_arr1[3] /= _full_prob;
            } else {
                for (int i = 0; i < 16; ++i) {
                    _nuc_arr1[i] /= _full_prob;
                }
                for (int i = 0; i < 16; ++i) {
                    _nuc_arr2[i] /= _full_prob;
                }
            }

            _pairs.reserve(_pair_map.size() + 40);
            for (auto it = _pair_map.begin(); it != _pair_map.end(); ++it) {
                _pairs.push_back(event_pair_t(it->first, it->second / _full_prob));
            }
            _pair_map.clear();
        }


        //
        // Forward-backward algorithm for VJ recombination receptors
        //

        // make a matrix chain with forward probabilities for VJ receptors
        void forward_vj(const MAAG &maag, event_ind_t j_ind) {
            this->fillZero(_forward_acc.get());

            // VJ probabilities for the fixed J
            for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_VAR_JOI_GEN_I); ++row_i) {
                (*_forward_acc)(VJ_VAR_JOI_GEN_I, 0, row_i, 0) = maag(VJ_VAR_JOI_GEN_I, 0, row_i, j_ind);
            }

            // V deletions
            for (event_ind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
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
        void backward_vj(const MAAG &maag, event_ind_t j_ind) {
            this->fillZero(_backward_acc.get());

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
            for (event_ind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_DEL_I); ++col_i) {
                    for (dim_t ins_col_i = 0; ins_col_i < maag.nodeColumns(VJ_VAR_JOI_INS_I); ++ins_col_i) {
                        (*_backward_acc)(VJ_VAR_DEL_I, v_ind, 0, col_i) +=
                                (*_backward_acc)(VJ_VAR_JOI_INS_I, 0, col_i, ins_col_i) * maag(VJ_VAR_JOI_INS_I, 0, col_i, ins_col_i);
                    }
                }
            }

            // V-J genes
            for (event_ind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                for (dim_t col_i = 0; col_i < maag.nodeColumns(VJ_VAR_DEL_I); ++col_i) {
                    (*_backward_acc)(VJ_VAR_JOI_GEN_I, 0, v_ind, 0) +=
                            (*_backward_acc)(VJ_VAR_DEL_I, v_ind, 0, col_i) * maag(VJ_VAR_DEL_I, v_ind, 0, col_i);
                }
            }

            for (event_ind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                _back_full_prob +=
                        (*_backward_acc)(VJ_VAR_JOI_GEN_I, 0, v_ind, 0) * maag(VJ_VAR_JOI_GEN_I, 0, v_ind, j_ind);
            }
        }


        // compute all forward-backward probabilities
        void forward_backward_vj(const MAAG &maag) {
            _forward_acc.reset(new ProbMMC());
            _forward_acc->resize(maag.chainSize());
            // VJ probabilities (for fixed J in future)
            _forward_acc->initNode(VJ_VAR_JOI_GEN_I, 1, maag.nodeRows(VJ_VAR_JOI_GEN_I), 1);
            // V deletions
            _forward_acc->initNode(VJ_VAR_DEL_I, maag.nodeSize(VJ_VAR_DEL_I), 1, maag.nodeColumns(VJ_VAR_DEL_I));
            // VJ insertions
            _forward_acc->initNode(VJ_VAR_JOI_INS_I, 1, maag.nodeRows(VJ_VAR_JOI_INS_I), maag.nodeColumns(VJ_VAR_JOI_INS_I));
            // J deletions (for fixed J in future)
            _forward_acc->initNode(VJ_JOI_DEL_I, 1, maag.nodeRows(VJ_JOI_DEL_I), 1);

            _backward_acc.reset(new ProbMMC());
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
            for (event_ind_t j_ind = 0; j_ind < maag.nJoi(); ++j_ind) {
                // compute forward and backward probabilities for a specific J gene
                this->forward_vj(maag, j_ind);
                this->backward_vj(maag, j_ind);

                // add fi * bi for this J to the accumulator
                for (dim_t row_i = 0; row_i < maag.nodeRows(VJ_VAR_JOI_GEN_I); ++row_i) {
                    this->pushEventPair(maag, VJ_VAR_JOI_GEN_I, 0, row_i, j_ind, 0, row_i, 0);
                }
                for (matrix_ind_t mat_i = 0; mat_i < maag.nVar(); ++mat_i) {
                    this->pushEventPairsWithErrors(maag, VJ_VAR_DEL_I, mat_i, mat_i, 0);
                }
                this->pushEventPairs(maag, VJ_VAR_JOI_INS_I, 0, 0);
                this->pushEventPairsWithErrors(maag, VJ_JOI_DEL_I, j_ind, 0, 1);

                this->inferInsertionNucleotides(maag, VJ_VAR_JOI_INS_I,
                                                0, maag.nodeColumns(VJ_VAR_DEL_I) - 1,
                                                maag.nodeColumns(VJ_VAR_DEL_I), maag.nodeColumns(VJ_VAR_DEL_I) + maag.nodeRows(VJ_JOI_DEL_I) - 1,
                                                _nuc_arr1);
            }

            this->vectorise_pair_map(maag);
        }


        //
        // Forward-backward algorithm for VDJ recombination receptors
        //

        // make a matrix chain with forward probabilities for VDJ receptors
        void forward_vdj(const MAAG &maag, event_ind_t d_ind, event_ind_t j_ind, bool recompute_d_gen_fi) {
            if (recompute_d_gen_fi) {
                // forward probabilities for (V prob -> V del -> VD ins) are fixed
                // for all pairs of J-D.
                // We have already stored in _forward_acc fi for this D, so we don't
                // need to recompute entire _forward_acc, we just need to recompute
                // J deletions and J genes fi.
                this->fillZero(_forward_acc.get(), VDJ_DIV_DEL_I);

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
                this->fillZero(_forward_acc.get(), VDJ_JOI_DEL_I);
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

            // update the full generation probability
            _full_prob += (*_forward_acc)(VDJ_JOI_DIV_GEN_I, 0, 0, 0);
        }


        // make a matrix chain with backward probabilities for VDJ receptors
        void backward_vdj(const MAAG &maag, event_ind_t d_ind, event_ind_t j_ind) {
            this->fillZero(_backward_acc.get());

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
            for (event_ind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
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
            _forward_acc.reset(new ProbMMC());
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

            _backward_acc.reset(new ProbMMC());
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
            this->fillZero(_forward_acc.get());
            // V genes and deletions
            for (event_ind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
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
                    for (event_ind_t v_ind = 0; v_ind < maag.nVar(); ++v_ind) {
                        (*_forward_acc)(VDJ_VAR_DIV_INS_I, 0, row_i, col_i) +=
                                (*_forward_acc)(VDJ_VAR_DEL_I, v_ind, 0, row_i) * maag(VDJ_VAR_DIV_INS_I, 0, row_i, col_i);
                    }
                }
            }

            // Compute fi * bi / Pgen for each event.
            bool recompute_d_gen_fi = true;
            for (event_ind_t d_ind = 0; d_ind < maag.nDiv(); ++d_ind) {
                recompute_d_gen_fi = true;
                for (event_ind_t j_ind = 0; j_ind < maag.nJoi(); ++j_ind) {
                    // compute forward and backward probabilities for a specific J gene
                    this->forward_vdj(maag, d_ind, j_ind, recompute_d_gen_fi);
                    this->backward_vdj(maag, d_ind, j_ind);
                    recompute_d_gen_fi = false;

                    // add fi * bi to the accumulator
                    for (matrix_ind_t mat_i = 0; mat_i < maag.nVar(); ++mat_i) {
                        this->pushEventPairs(maag, VDJ_VAR_GEN_I, mat_i, mat_i);
                        this->pushEventPairsWithErrors(maag, VDJ_VAR_DEL_I, mat_i, mat_i, 0);
                    }
                    this->pushEventPairs(maag, VDJ_VAR_DIV_INS_I, 0, 0);
                    this->pushEventPairsWithErrors(maag, VDJ_DIV_DEL_I, d_ind, 0, 1);
                    this->pushEventPairs(maag, VDJ_DIV_JOI_INS_I, 0, 0);
                    this->pushEventPairsWithErrors(maag, VDJ_JOI_DEL_I, j_ind, 0, 2);
                    this->pushEventPair(maag, VDJ_JOI_DIV_GEN_I, 0, j_ind, d_ind, 0, 0, 0);

                    seq_len_t v_vertices = maag.nodeColumns(VDJ_VAR_DEL_I),
                            d3_vertices = maag.nodeRows(VDJ_DIV_DEL_I),
                            d5_vertices = maag.nodeColumns(VDJ_DIV_DEL_I),
                            j_vertices = maag.nodeRows(VDJ_JOI_DEL_I);

                    this->inferInsertionNucleotides(maag, VDJ_VAR_DIV_INS_I,
                                                    0, v_vertices - 1,
                                                    v_vertices, v_vertices + d3_vertices - 1,
                                                    _nuc_arr1, false);


                    this->inferInsertionNucleotides(maag, VDJ_DIV_JOI_INS_I,
                                                    v_vertices + d3_vertices, v_vertices + d3_vertices + d5_vertices - 1,
                                                    v_vertices + d3_vertices + d5_vertices, v_vertices + d3_vertices + d5_vertices + j_vertices - 1,
                                                    _nuc_arr2, true);
                }
            }

            this->vectorise_pair_map(maag);
        }

    };


    class ForwardBackwardAlgorithm {
    public:


        ForwardBackwardAlgorithm() { }


        virtual ~ForwardBackwardAlgorithm() { }


        virtual void infer(const MAAG &maag) = 0;


        event_pair_t nextEvent() {
            if (_status) {
                event_pair_t res = _pairs[_pairs_i];
//                cout << res.second << " -> ";
//                res.second = res.second / _full_prob;
//                cout << res.second << endl;
                ++_pairs_i;
                return res;
            }

            _status = false;
            return event_pair_t(0, 0);
        }


        bool is_empty() {
//            return _status && _pairs_i != _pairs.size();
            if (_status) {
                if (_pairs_i != _pairs.size()) {
                    return false;
                }
            }
            return true;
        }


        bool status() const { return _status; }


        prob_t fullProbability() const { return _full_prob; }


        prob_t bfullProbability() const { return _back_full_prob; }


        const vector<event_pair_t>& event_pairs() const { return _pairs; }


        prob_t* insertion_probs() const;


    protected:
        pProbMMC _forward_acc, _backward_acc;  /** Temporary MMC for storing forward and backward probabilities correspondingly. */
        bool _status;
        vector<event_pair_t> _pairs;
        size_t _pairs_i;
        prob_t _full_prob, _back_full_prob;

    };


    class VJForwardBackward : public ForwardBackwardAlgorithm {

    };


    class VDJForwardBackward : public ForwardBackwardAlgorithm {

    };


    class VD2JForwardBackward : public ForwardBackwardAlgorithm {

    };


}

#endif //YMIR_MAAGFORWARDBACKWARDALGORITHM_H
