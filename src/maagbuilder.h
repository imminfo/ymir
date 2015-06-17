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

#ifndef _MAAGBUILDER_H_
#define _MAAGBUILDER_H_

#define DEFAULT_SEQ_POSES_RESERVE 300

#define VARIABLE_GENES_MATRIX_INDEX 0
#define VARIABLE_DELETIONS_MATRIX_INDEX 1
#define JOINING_DELETIONS_VJ_MATRIX_INDEX 3
#define JOINING_GENES_VDJ_MATRIX_INDEX 6
#define JOINING_DELETIONS_VDJ_MATRIX_INDEX 5
#define DIVERSITY_GENES_MATRIX_INDEX 3
#define VarJoi_INSERTIONS_MATRIX_INDEX 2
#define VarDiv_INSERTIONS_MATRIX_INDEX 2
#define DivJoi_INSERTIONS_MATRIX_INDEX 4


#include "maag.h"
#include "genesegment.h"
#include "modelparametervector.h"
#include "markovchain.h"


namespace ymir {

    class MAAGBuilder;


    /**
    * \class MAAGBuilder
    */
    class MAAGBuilder : protected MAAG {

    public:


        MAAGBuilder(const ModelParameterVector &param_vec, const VDJRecombinationGenes &genes) : MAAG() {
            _param_vec = new ModelParameterVector(param_vec);
            _genes = new VDJRecombinationGenes(genes);
        }


        virtual ~MAAGBuilder() {
            if (_param_vec) { delete _param_vec; }
            if (_genes) { delete _genes; }
        }


        /**
         * \brief Build MAAGs from the given clonotypes.
         *
         * \param clonotype Clonotype from which build the MAAG.
         * \param cloneset Set of clonotypes from which build the repertoires of MAAGs.
         * \param full_build If true than make MAAG with stored event indices.
         * ???? \param aminoacid If true than build MAAGs from the aminoacid sequences.
         */
        ///@{
        MAAG build(const Clonotype &clonotype, bool full_build = false) const {
            ProbMMC probs;
            EventIndMMC events;
            vector<seq_len_t> seq_poses;
            seq_poses.reserve(DEFAULT_SEQ_POSES_RESERVE);

            if (clonotype.is_vj()) {
                probs.resize(VJ_CHAIN_SIZE);
                if (full_build) { events.resize(VJ_CHAIN_SIZE); }
            } else {
                probs.resize(VDJ_CHAIN_SIZE);
                if (full_build) { events.resize(VDJ_CHAIN_SIZE); }
            }

            this->buildVariable(clonotype, probs, events, seq_poses, full_build);
            this->buildJoining(clonotype, probs, events, seq_poses, full_build);
            if (clonotype.is_vj()) {
                this->buildVJinsertions(clonotype, probs, events, seq_poses, full_build);
            } else {
                this->buildDiversity(clonotype, probs, events, seq_poses, full_build);
                this->buildVDinsertions(clonotype, probs, events, seq_poses, full_build);
                this->buildDJinsertions(clonotype, probs, events, seq_poses, full_build);
            }

            probs.finish();
            events.finish();
            if (full_build) {

                // TODO:
                // if DJ insertions matrix is consists only of zeros for gene Dn, then
                // remove VD insertions, D deletions and DJ insertions for the Dn gene.
                //

                seq_len_t *seq_poses_arr = new seq_len_t[seq_poses.size()];
                copy(seq_poses.begin(), seq_poses.end(), seq_poses_arr);
                return MAAG(probs, events, clonotype.sequence(), seq_poses_arr, seq_poses.size());
            } else {
                return MAAG(probs);
            }
        }


        MAAG* buildPtr(const Clonotype &clonotype, bool full_build = false) const {
            ProbMMC probs;
            EventIndMMC events;
            vector<seq_len_t> seq_poses;
            seq_poses.reserve(DEFAULT_SEQ_POSES_RESERVE);

            if (clonotype.is_vj()) {
                probs.resize(VJ_CHAIN_SIZE);
                if (full_build) { events.resize(VJ_CHAIN_SIZE); }
            } else {
                probs.resize(VDJ_CHAIN_SIZE);
                if (full_build) { events.resize(VDJ_CHAIN_SIZE); }
            }

            this->buildVariable(clonotype, probs, events, seq_poses, full_build);
            this->buildJoining(clonotype, probs, events, seq_poses, full_build);
            if (clonotype.is_vj()) {
                this->buildVJinsertions(clonotype, probs, events, seq_poses, full_build);
            } else {
                this->buildDiversity(clonotype, probs, events, seq_poses, full_build);
                this->buildVDinsertions(clonotype, probs, events, seq_poses, full_build);
                this->buildDJinsertions(clonotype, probs, events, seq_poses, full_build);
            }

            probs.finish();
            events.finish();
            if (full_build) {
                seq_len_t *seq_poses_arr = new seq_len_t[seq_poses.size()];
                copy(seq_poses.begin(), seq_poses.end(), seq_poses_arr);
                return new MAAG(probs, events, clonotype.sequence(), seq_poses_arr, seq_poses.size());
            } else {
                return new MAAG(probs);
            }
        }


        MAAGRepertoire build(const ClonesetView &cloneset, bool full_build = false) const {
            MAAGRepertoire res;
//            res.reserve(cloneset.size());
            res.resize(cloneset.size());
            MAAG tmp;
            for (size_t i = 0; i < cloneset.size(); ++i) {
//                res.push_back(MAAG(&(this->build(cloneset[i], full_build))));

                res[i] = this->build(cloneset[i], full_build);

//                tmp = this->build(cloneset[i], full_build);
//                res[i].swap_maag(this->build(cloneset[i], full_build));

                if ((i+1) % 50000 == 0) {
                    cout << "Built " << (int) (i+1) << " graphs." << endl;
                }
            }
            return res;
        }
        ///@}


        prob_t buildAndCompute(const Clonotype &clonotype, bool aminoacid = false, MAAG_COMPUTE_PROB_ACTION action = SUM_PROBABILITY) const {
            return this->build(clonotype, false).fullProbability(action);
        }

        vector<prob_t> buildAndCompute(const ClonesetView &cloneset, bool aminoacid = false, MAAG_COMPUTE_PROB_ACTION action = SUM_PROBABILITY) const {
            vector<prob_t> res;
            res.reserve(cloneset.size());
            for (size_t i = 0; i < cloneset.size(); ++i) {
                res.push_back(buildAndCompute(cloneset[i], aminoacid, action));
                if ((i+1) % 50000 == 0) {
                    cout << "Computed " << (int) (i+1) << " assembling probabilities." << endl;
                }
            }
            return res;
        }


        /**
         * \brief Replace event probabilities in the given MAAGs if they have stored event indices.
         */
        ///@{
        void updateEventProbabilities(MAAG *maag) const {
            if (maag->has_events()) {
                for (int node_i = 0; node_i < maag->chainSize(); ++node_i) {
                    // either rebuild all insertions
                    if (maag->is_vj() && node_i == VarJoi_INSERTIONS_MATRIX_INDEX) {
                        // ???
                    } else if (maag->is_vdj() && node_i == VarDiv_INSERTIONS_MATRIX_INDEX) {
                        // ???
                    } else if (maag->is_vdj() && node_i == DivJoi_INSERTIONS_MATRIX_INDEX) {
                        // ???
                    } else {
                        // or just replace all event probabilities with the new ones
                        for (int mat_i = 0; mat_i < maag->nodeSize(node_i); ++mat_i) {
                            for (int row_i = 0; row_i < maag->nodeRows(node_i); ++row_i) {
                                for (int col_i = 0; col_i < maag->nodeColumns(node_i); ++col_i) {
                                    (*maag)(node_i, mat_i, row_i, col_i) =
                                            (*_param_vec)[maag->event_index(node_i, mat_i, row_i, col_i)];
                                }
                            }
                        }
                    }
                }
            }
        }

        void updateEventProbabilities(MAAGRepertoire *repertoire) const {
            for (size_t i = 0; i < repertoire->size(); ++i) {
                this->updateEventProbabilities(repertoire->data() + i);
            }
        }
        ///@}


    protected:

        ModelParameterVector *_param_vec;  // or just copy it?
        VDJRecombinationGenes *_genes; // copy this too?


        /**
         * \brief Private default constructor.
         */
        MAAGBuilder() : _param_vec(nullptr), _genes(nullptr) {}


        /**
        * \brief Build probability and events matrices for Variable gene segments.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param full_build Boolean if build should be full.
        */
        void buildVariable(const Clonotype &clonotype,
                           ProbMMC &probs,
                           EventIndMMC &events,
                           vector<seq_len_t> &seq_poses,
                           bool full_build) const
        {
            // find max V alignment
            seq_len_t len = 0;
            segindex_t v_num = clonotype.nVar(), j_num = clonotype.nJoi();
            for (int v_index = 0; v_index < v_num; ++v_index) {
                if (clonotype.getVend(v_index) > len) {
                    len = clonotype.getVend(v_index);
                }
            }

            // compute V deletions
            seq_len_t v_len = 0;
            segindex_t v_gene = 0;
            seq_len_t v_end = 0;

            probs.initNode(VARIABLE_DELETIONS_MATRIX_INDEX, v_num, 1, len + 1);
            if (full_build) {
                events.initNode(VARIABLE_DELETIONS_MATRIX_INDEX, v_num, 1, len + 1);
            }

            if (clonotype.is_vj()) {
                probs.initNode(VARIABLE_GENES_MATRIX_INDEX, 1, v_num, j_num);
                if (full_build) {
                    events.initNode(VARIABLE_GENES_MATRIX_INDEX, 1, v_num, j_num);
                }
            } else {
                probs.initNode(VARIABLE_GENES_MATRIX_INDEX, v_num, 1, 1);
                if (full_build) {
                    events.initNode(VARIABLE_GENES_MATRIX_INDEX, v_num, 1, 1);
                }
            }

            EVENT_CLASS V_DEL = clonotype.is_vj() ? VJ_VAR_DEL : VDJ_VAR_DEL;
            for (segindex_t v_index = 0; v_index < v_num; ++v_index) {
                v_gene = clonotype.getVar(v_index);
                v_len = _genes->V()[v_gene].sequence.size();
                v_end = clonotype.getVend(v_index);

                if (clonotype.is_vj()) {
                    // probability of choosing this V gene segment
                    for (segindex_t j_index = 0; j_index < j_num; ++j_index) {
                        probs(VARIABLE_GENES_MATRIX_INDEX, 0, v_index, j_index)
                                = _param_vec->event_prob(VJ_VAR_JOI_GEN, 0, v_gene - 1, clonotype.getJoi(j_index) - 1);
                    }
                } else {
                    probs(VARIABLE_GENES_MATRIX_INDEX, v_index, 0, 0) = _param_vec->event_prob(VDJ_VAR_GEN, 0, v_gene - 1); // probability of choosing this V gene segment
                }

                for (seq_len_t i = 0; i < len + 1; ++i) {
                    if (v_len - i >= 0 && i <= v_end) {
                        probs(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = _param_vec->event_prob(V_DEL, v_gene - 1, v_len - i); // probability of deletions
                    } else {
                        probs(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = 0; // if exceeds length of V gene segment
                    }
                    // probs(1, v_index, 0, i) = (v_len - i >= 0) ? _param_vec->prob_V_del(v_gene, v_len - i) : 0;
                    // OPTIMISATION: first find where zeros are start and then just fill this part of the vector with zeros without any "ifs"
                }

                if (full_build) {
                    if (clonotype.is_vj()) {
                        // probability of choosing this V gene segment
                        for (segindex_t j_index = 0; j_index < j_num; ++j_index) {
                            events(VARIABLE_GENES_MATRIX_INDEX, 0, v_index, j_index)
                                    = _param_vec->event_index(VJ_VAR_JOI_GEN, 0, v_gene - 1, clonotype.getJoi(j_index) - 1);
                        }
                    } else {
                        events(VARIABLE_GENES_MATRIX_INDEX, v_index, 0, 0) = _param_vec->event_index(VDJ_VAR_GEN, 0, v_gene - 1);
                    }

                    for (seq_len_t i = 0; i < len + 1; ++i) {
                        if (v_len - i >= 0 && i <= v_end) {
                            events(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = _param_vec->event_index(V_DEL, v_gene - 1, v_len - i);
                        } else {
                            events(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = 0;
                        }
                    }
                }
            }

            for (seq_len_t i = 0; i < len + 1; ++i) {
                seq_poses.push_back(i);
            }
        }


        /**
        * \brief Build probability and events matrices for Joining gene segments.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param full_build Boolean if build should be full.
        */
        void buildJoining(const Clonotype &clonotype,
                          ProbMMC &probs,
                          EventIndMMC &events,
                          vector<seq_len_t> &seq_poses,
                          bool full_build) const
        {
            int J_index_dels = JOINING_DELETIONS_VJ_MATRIX_INDEX,
                    J_index_genes = JOINING_GENES_VDJ_MATRIX_INDEX;
            if (clonotype.is_vdj()) {
                J_index_dels = JOINING_DELETIONS_VDJ_MATRIX_INDEX;
            }

            // find max J alignment
            segindex_t j_num = clonotype.nJoi();
            seq_len_t len = clonotype.sequence().size();
            for (int j_index = 0; j_index < j_num; ++j_index) {
                if (clonotype.getJstart(j_index) < len) {
                    len = clonotype.getJstart(j_index);
                }
            }
            len = clonotype.sequence().size() - len + 1;


            // add J deletions nodes
            probs.initNode(J_index_dels, j_num, len + 1, 1);
            if (full_build) {
                events.initNode(J_index_dels, j_num, len + 1, 1);
            }

            // add J or J-D gene nodes
            if (clonotype.is_vdj()) {
                probs.initNode(J_index_genes, 1, j_num, clonotype.nDiv());
                if (full_build) {
                    events.initNode(J_index_genes, 1, j_num, clonotype.nDiv());
                }
            }

            // compute J deletions
            seq_len_t j_len = 0;
            segindex_t j_gene = 0;
            seq_len_t j_start = 0;

            EVENT_CLASS J_DEL = clonotype.is_vj() ? VJ_JOI_DEL : VDJ_JOI_DEL;
            for (segindex_t j_index = 0; j_index < j_num; ++j_index) {
                j_gene = clonotype.getJoi(j_index);
                j_len = _genes->J()[j_gene].sequence.size();
                j_start = clonotype.getJstart(j_index);

                if (clonotype.is_vdj()) {
                    for (segindex_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
                        probs(J_index_genes, 0, j_index, d_index)
                                = _param_vec->event_prob(VDJ_JOI_DIV_GEN, 0, j_gene - 1, clonotype.getDiv(d_index) - 1); // probability of choosing this J gene segment with other D genes
                    }
                }

                for (seq_len_t i = 0; i < len + 1; ++i) {
                    if (j_len - len + i >= 0 && len - i <= clonotype.sequence().size() - j_start + 1) {
                        probs(J_index_dels, j_index, i, 0) = _param_vec->event_prob(J_DEL, j_gene - 1, j_len - len + i); // probability of deletions
                    } else {
                        probs(J_index_dels, j_index, i, 0) = 0; // if exceeds length of J gene segment
                    }

                }

                if (full_build) {
                    if (clonotype.is_vdj()) {
                        for (segindex_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
                            events(J_index_genes, 0, j_index, d_index)
                                    = _param_vec->event_index(VDJ_JOI_DIV_GEN, 0, j_gene - 1, clonotype.getDiv(d_index) - 1); // probability of choosing this J gene segment with other D genes
                        }
                    }

                    for (seq_len_t i = 0; i < len + 1; ++i) {
                        if (j_len - len + i >= 0 && len - i <= clonotype.sequence().size() - j_start + 1) {
                            events(J_index_dels, j_index, i, 0) = _param_vec->event_index(J_DEL, j_gene - 1, j_len - len + i);
                        } else {
                            events(J_index_dels, j_index, i, 0) = 0;
                        }
                    }
                }
            }

            for (seq_len_t i = clonotype.sequence().size() - len + 1; i < clonotype.sequence().size() + 1; ++i) {
                seq_poses.push_back(i);
            }
        }


        /**
        * \brief Build probability and events matrices for Diversity gene segments.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param full_build Boolean if build should be full.
        */
        void buildDiversity(const Clonotype &clonotype,
                            ProbMMC &probs,
                            EventIndMMC &events,
                            vector<seq_len_t> &seq_poses,
                            bool full_build) const
        {
            d_alignment_t d_alignment;

            // vector seq_start -> 0 means no such index in the matrix, 1 otherwise.
            seq_len_t seq_arr_size = clonotype.sequence().size() + 1;
            seq_len_t *seq_row = new seq_len_t[seq_arr_size];
            std::fill(seq_row, seq_row + seq_arr_size, 0);
            // vector seq_end -> 0 means no such index in the matrix, 1 otherwise.
            seq_len_t *seq_col = new seq_len_t[seq_arr_size];
            std::fill(seq_col, seq_col + seq_arr_size, 0);

            for (segindex_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
                seq_len_t min_D_len = _param_vec->D_min_len(clonotype.getDiv(d_index));

                for (segindex_t j = 0; j < clonotype.nDalignments(d_index); ++j) {
                    d_alignment = clonotype.getDalignment(d_index, j);

                    // yes-yes, I know that it could be done more efficiently. But I don't want to.
                    for (seq_len_t i = d_alignment.seqstart; i <= d_alignment.seqend - min_D_len + 1; ++i) {
                        seq_row[i] = 1;
                    }

                    for (seq_len_t i = d_alignment.seqstart + min_D_len - (seq_len_t) 1; i <= d_alignment.seqend; ++i) {
                        seq_col[i] = 1;
                    }
                }
            }

            // make new vector seq_start -> 1based index in rows of Ddel matrices
            seq_len_t seq_ind = 1;
            for (seq_len_t i = 0; i < seq_arr_size; ++i) {
                if (seq_row[i]) {
                    seq_row[i] = seq_ind;
                    ++seq_ind;
                }
            }

            // make new vector seq_end -> 1based index in columns of Ddel matrices
            seq_ind = 1;
            for (seq_len_t i = 0; i < seq_arr_size; ++i) {
                if (seq_col[i]) {
                    seq_col[i] = seq_ind;
                    ++seq_ind;
                }
            }

            // find indices of D alignments and use only them to reduce memory usage.
            seq_len_t last_max_seq_start = 0, last_max_seq_end = 0, seq_row_ind = 0, seq_col_ind = 0;

            seq_len_t seq_row_nonzeros = 0, seq_col_nonzeros = 0;
            for (seq_len_t i = 0; i < seq_arr_size; ++i) {
                seq_row_nonzeros += seq_row[i] != 0;
                seq_col_nonzeros += seq_col[i] != 0;
            }

            probs.initNode(DIVERSITY_GENES_MATRIX_INDEX, clonotype.nDiv(), seq_row_nonzeros, seq_col_nonzeros);
            if (full_build) {
                events.initNode(DIVERSITY_GENES_MATRIX_INDEX, clonotype.nDiv(), seq_row_nonzeros, seq_col_nonzeros);
            }


            segindex_t d_index = 0, d_gene = 0;
            seq_len_t min_D_len = 0, d_len = 0;

            for (segindex_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
                d_gene = clonotype.getDiv(d_index);
                d_len = _genes->D()[d_gene].sequence.size();
                min_D_len = _param_vec->D_min_len(d_gene);

                // for each aligned Div segment get all possible smaller alignments and add them to the matrix.
                for (segindex_t j = 0; j < clonotype.nDalignments(d_index); ++j) {
                    d_alignment = clonotype.getDalignment(d_index, j);

                    for (seq_len_t left_pos = d_alignment.seqstart; left_pos <= d_alignment.seqend - min_D_len + 1; ++left_pos) {
                        for (seq_len_t right_pos = left_pos + min_D_len - 1; right_pos <= d_alignment.seqend; ++right_pos) {
                            probs(DIVERSITY_GENES_MATRIX_INDEX, d_index, seq_row[left_pos] - 1, seq_col[right_pos] - 1)
                                    = _param_vec->event_prob(VDJ_DIV_DEL,
                                                             d_gene - 1,
                                                             d_alignment.Dstart + left_pos - d_alignment.seqstart - 1,
                                                             d_len - (d_alignment.Dend - (d_alignment.seqend - right_pos)));
                            if (full_build) {
                                events(DIVERSITY_GENES_MATRIX_INDEX, d_index, seq_row[left_pos] - 1, seq_col[right_pos] - 1)
                                        = _param_vec->event_index(VDJ_DIV_DEL,
                                                                  d_gene - 1,
                                                                  d_alignment.Dstart + left_pos - d_alignment.seqstart - 1,
                                                                  d_len - (d_alignment.Dend - (d_alignment.seqend - right_pos)));
                            }
                        }
                    }
                }
            }

            // insert D3 and D5 positions
            vector<seq_len_t> D35_poses;
            D35_poses.reserve(seq_row_nonzeros + seq_col_nonzeros + 2);
            for (seq_len_t i = 1; i < seq_arr_size; ++i) { if (seq_row[i]) { D35_poses.push_back(i); } }
            for (seq_len_t i = 1; i < seq_arr_size; ++i) { if (seq_col[i]) { D35_poses.push_back(i); } }

            delete [] seq_row;
            delete [] seq_col;

            // Note! insert diversity gene seq poses BEFORE joining gene seq poses
            seq_poses.reserve(seq_poses.size() + D35_poses.size() + 2);  // +2 -> just in case (:
            seq_poses.insert(seq_poses.begin() + probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX), D35_poses.begin(), D35_poses.end());
        }


        /**
        * \brief Build probability and events matrices for Variable-Joining gene segments insertions.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param full_build Boolean if build should be full.
        */
        void buildVJinsertions(const Clonotype &clonotype,
                               ProbMMC &probs,
                               EventIndMMC &events,
                               vector<seq_len_t> &seq_poses,
                               bool full_build) const
        {
            MarkovChain mc(_param_vec->get_iterator(_param_vec->event_index(VJ_VAR_JOI_INS_NUC, 0, 0)));

            seq_len_t v_vertices = probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                    j_vertices = probs.nodeRows(JOINING_DELETIONS_VJ_MATRIX_INDEX);

            probs.initNode(VarJoi_INSERTIONS_MATRIX_INDEX, 1, v_vertices, j_vertices);

            if (full_build) {
                events.initNode(VarJoi_INSERTIONS_MATRIX_INDEX, 1, v_vertices, j_vertices);
            }

            this->buildInsertions(clonotype,
                                  probs,
                                  events,
                                  seq_poses,
                                  VarJoi_INSERTIONS_MATRIX_INDEX,
                                  _param_vec->event_index(VJ_VAR_JOI_INS_LEN, 0, 0),
                                  _param_vec->max_VJ_ins_len(),
                                  full_build,
                                  0,
                                  v_vertices - 1,
                                  v_vertices,
                                  v_vertices + j_vertices - 1,
                                  mc);
        }


        /**
        * \brief Build probability and events matrices for Variable-Diversity gene segments insertions.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param full_build Boolean if build should be full.
        */
        void buildVDinsertions(const Clonotype &clonotype,
                               ProbMMC &probs,
                               EventIndMMC &events,
                               vector<seq_len_t> &seq_poses,
                               bool full_build) const
        {
            MarkovChain mc(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC, 0, 0)));

            seq_len_t v_vertices = probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                    d3_vertices = probs.nodeRows(DIVERSITY_GENES_MATRIX_INDEX);

            probs.initNode(VarDiv_INSERTIONS_MATRIX_INDEX, 1, v_vertices, d3_vertices);

            if (full_build) {
                events.initNode(VarDiv_INSERTIONS_MATRIX_INDEX, 1, v_vertices, d3_vertices);
            }

            this->buildInsertions(clonotype,
                                  probs,
                                  events,
                                  seq_poses,
                                  VarDiv_INSERTIONS_MATRIX_INDEX,
                                  _param_vec->event_index(VDJ_VAR_DIV_INS_LEN, 0, 0),
                                  _param_vec->max_VD_ins_len(),
                                  full_build,
                                  0,
                                  v_vertices - 1,
                                  v_vertices,
                                  v_vertices + d3_vertices - 1,
                                  mc);
        }


        /**
        * \brief Build probability and events matrices for Diversity-Joining gene segments insertions.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param full_build Boolean if build should be full.
        */
        void buildDJinsertions(const Clonotype &clonotype,
                               ProbMMC &probs,
                               EventIndMMC &events,
                               vector<seq_len_t> &seq_poses,
                               bool full_build) const
        {
            MarkovChain mc(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC, 0, 0)));

            seq_len_t v_vertices = probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                    d3_vertices = probs.nodeRows(DIVERSITY_GENES_MATRIX_INDEX),
                    d5_vertices = probs.nodeColumns(DIVERSITY_GENES_MATRIX_INDEX),
                    j_vertices = probs.nodeRows(JOINING_DELETIONS_VDJ_MATRIX_INDEX);

            probs.initNode(DivJoi_INSERTIONS_MATRIX_INDEX, 1, d5_vertices, j_vertices);

            if (full_build) {
                events.initNode(DivJoi_INSERTIONS_MATRIX_INDEX, 1, d5_vertices, j_vertices);
            }

            this->buildInsertions(clonotype,
                                  probs,
                                  events,
                                  seq_poses,
                                  DivJoi_INSERTIONS_MATRIX_INDEX,
                                  _param_vec->event_index(VDJ_DIV_JOI_INS_LEN, 0, 0),
                                  _param_vec->max_DJ_ins_len(),
                                  full_build,
                                  v_vertices + d3_vertices,
                                  v_vertices + d3_vertices + d5_vertices - 1,
                                  v_vertices + d3_vertices + d5_vertices,
                                  v_vertices + d3_vertices + d5_vertices + j_vertices - 1,
                                  mc);
        }


        /**
         * \brief General function for building insertions.
         *
         * \param clonotype Clonotype that used for building the graph.
         * \param probs Multi-Matrix Chain with event probabilities.
         * \param events Multi-Matrix Chain with event indices.
         * \param seq_poses Vector of positions.
         * \param ins_node_index Node of the event / prob matrix with insertions.
         * \param null_insertion Event index of insertions of length zero.
         * \param max_size Maximum size of the length of insertions.
         * \param full_build Boolean if build should be full.
         * \param left_vertices_start Starting index in seq_poses for the vertices in the left matrix.
         * \param left_vertices_end Ending index in seq_poses for the vertices in the left matrix.
         * \param right_vertices_start Starting index in seq_poses for the vertices in the right matrix.
         * \param right_vertices_end Ending index in seq_poses for the vertices in the right matrix.
         * \param mc MarkovChain that use for generation of N nucleotides.
        */
        void buildInsertions(const Clonotype &clonotype,
                             ProbMMC &probs,
                             EventIndMMC &events,
                             vector<seq_len_t> &seq_poses,
                             seq_len_t ins_node_index,
                             eventind_t null_insertion,
                             seq_len_t max_size,
                             bool full_build,
                             seq_len_t left_vertices_start,
                             seq_len_t left_vertices_end,
                             seq_len_t right_vertices_start,
                             seq_len_t right_vertices_end,
                             const MarkovChain& mc) const
        {
            int insertion_len;
            bool good_insertion;

            for (size_t left_vertex_i = left_vertices_start; left_vertex_i <= left_vertices_end; ++left_vertex_i) {
                for (size_t right_vertex_i = right_vertices_start; right_vertex_i <= right_vertices_end; ++right_vertex_i) {
                    insertion_len = seq_poses[right_vertex_i] - seq_poses[left_vertex_i] - 1;
                    good_insertion = (insertion_len >= 0) && (insertion_len <= max_size);
                    if (good_insertion) {
                        if (seq_poses[left_vertex_i] == 0) {
                            // seq_iterator(0) and multiplication by .25 here because we don't know the last nucleotide
                            // of the aligned gene segment. And, honestly, it doesn't even matter and won't dramatically affect the result.
                            probs(ins_node_index, 0, left_vertex_i - left_vertices_start, right_vertex_i - right_vertices_start)
                                    = mc.nucProbability(clonotype.seq_iterator(0), insertion_len)
                                      * .25
                                      * (*_param_vec)[null_insertion + insertion_len];
                        } else {
                            probs(ins_node_index, 0, left_vertex_i - left_vertices_start, right_vertex_i - right_vertices_start)
                                    = mc.nucProbability(clonotype.seq_iterator(seq_poses[left_vertex_i] - 1), insertion_len)
                                      * (*_param_vec)[null_insertion + insertion_len];
                        }

                        if (full_build) {
                            events(ins_node_index, 0, left_vertex_i - left_vertices_start, right_vertex_i - right_vertices_start)
                                    = null_insertion + insertion_len;
                        }
                    } else {
                        probs(ins_node_index, 0, left_vertex_i - left_vertices_start, right_vertex_i - right_vertices_start) = 0;
                        if (full_build) {
                            events(ins_node_index, 0, left_vertex_i - left_vertices_start, right_vertex_i - right_vertices_start) = 0;
                        }
                    }
                }
            }
        }

    };
}

#endif //_MAAGBUILDER_H_
