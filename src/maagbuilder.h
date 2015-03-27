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

#define VJ_CHAIN_SIZE 5
#define VDJ_CHAIN_SIZE 7

#define VARIABLE_GENES_MATRIX_INDEX 0
#define VARIABLE_DELETIONS_MATRIX_INDEX 1
#define JOINING_GENES_VJ_MATRIX_INDEX 3
#define JOINING_DELETIONS_VJ_MATRIX_INDEX 4
#define JOINING_GENES_VDJ_MATRIX_INDEX 5
#define JOINING_DELETIONS_VDJ_MATRIX_INDEX 6
#define DIVERSITY_GENES_MATRIX_INDEX 3
#define VarJoi_INSERTIONS_MATRIX_INDEX 2
#define VarDiv_INSERTIONS_MATRIX_INDEX 2
#define DivJoi_INSERTIONS_MATRIX_INDEX 4


#include "maag.h"


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


        MAAG build(/*const*/ Clonotype &clonotype, bool full_build = false) const {
            ProbMMC probs;
            EventIndMMC events;
            vector<seq_len_t> seq_poses;
            seq_poses.reserve(DEFAULT_SEQ_POSES_RESERVE);


            if (!clonotype.is_vdj()) {
                probs.resize(VJ_CHAIN_SIZE);
                if (full_build) { events.resize(VJ_CHAIN_SIZE); }
            } else {
                probs.resize(VDJ_CHAIN_SIZE);
                if (full_build) { events.resize(VDJ_CHAIN_SIZE); }
            }

            this->buildVarible(clonotype, probs, events, seq_poses, full_build);
            this->buildJoining(clonotype, probs, events, seq_poses, full_build);
            if (!clonotype.is_vdj()) {
                this->buildVJinsertions(clonotype, probs, events, seq_poses, full_build);
            } else {
                this->buildDiversity(clonotype, probs, events, seq_poses, full_build);
                this->buildVDinsertions(clonotype, probs, events, seq_poses, full_build);
                this->buildDJinsertions(clonotype, probs, events, seq_poses, full_build);
            }

            seq_len_t *seq_poses_arr = new seq_len_t[seq_poses.size()];
            copy(seq_poses.begin(), seq_poses.end(), seq_poses_arr);
//            for (int i = 0; i < seq_poses.size(); ++i) { seq_poses_arr[i] = seq_poses[i]; }

            if (full_build) {
                return MAAG(probs, events, clonotype.sequence(), seq_poses_arr, seq_poses.size());
            } else {
                return MAAG(probs);
            }
        }


//        MAAGRepertoire build(const Cloneset &cloneset, bool full_build = false) {
//            // ???
//        }


    protected:

        ModelParameterVector *_param_vec;  // or just copy it?
        VDJRecombinationGenes *_genes; // copy this too?
        // MarkovChain *mc_alpha, *mc_beta;


        MAAGBuilder() {}


        /**
        * \brief Build probability and events matrices for Variable gene segments.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param full_build Boolean if build should be full.
        */
        void buildVarible(const Clonotype &clonotype,
                ProbMMC &probs,
                EventIndMMC &events,
                vector<seq_len_t> &seq_poses,
                bool full_build) const
        {
            // find max V alignment
            seq_len_t len = 0;
            segindex_t v_num = clonotype.nV();
            for (int v_index = 0; v_index < v_num; ++v_index) {
                if (clonotype.getVend(v_index) > len) {
                    len = clonotype.getVend(v_index);
                }
            }

            // compute V deletions
            seq_len_t v_len = 0;
            segindex_t v_gene = 0;
            seq_len_t v_end = 0;
            probs.initNode(VARIABLE_GENES_MATRIX_INDEX, v_num, 1, 1);
            probs.initNode(VARIABLE_DELETIONS_MATRIX_INDEX, v_num, 1, len + 1);
            if (full_build) {
                events.initNode(VARIABLE_GENES_MATRIX_INDEX, v_num, 1, 1);
                events.initNode(VARIABLE_DELETIONS_MATRIX_INDEX, v_num, 1, len + 1);
            }
            for (segindex_t v_index = 0; v_index < v_num; ++v_index) {
                v_gene = clonotype.getV(v_index);
                v_len = _genes->V()[v_gene].sequence.size();
                v_end = clonotype.getVend(v_index);

                probs(VARIABLE_GENES_MATRIX_INDEX, v_index, 0, 0) = _param_vec->prob_V_gene(v_gene); // probability of choosing this V gene segment
                for (seq_len_t i = 0; i < len + 1; ++i) {
                    if (v_len - i >= 0 && v_len - i <= v_end) {
                        probs(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = _param_vec->prob_V_del(v_gene, v_len - i); // probability of deletions
                    } else {
                        probs(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = 0; // if exceeds length of V gene segment
                    }
                    // probs(1, v_index, 0, i) = (v_len - i >= 0) ? _param_vec->prob_V_del(v_gene, v_len - i) : 0;
                    // OPTIMISATION: first find where zeros are start and then just fill this part of the vector with zeros without any "ifs"
                }

                if (full_build) {
                    events(VARIABLE_GENES_MATRIX_INDEX, v_index, 0, 0) = _param_vec->index_V_gene(v_gene);
                    for (seq_len_t i = 0; i < len + 1; ++i) {
                        if (v_len - i >= 0 && v_len - i <= v_end) {
                            events(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = _param_vec->index_V_del(v_gene, v_len - i);
                        } else {
                            events(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = 0;
                        }
                    }

                    for (seq_len_t i = 0; i < len + 1; ++i) {
                        seq_poses.push_back(i);
                    }
                }
            }

        }


        /**
        * \brief Build probability and events matrices for Variable gene segments.
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
            int J_index_dels = probs.chainSize() - 2, J_index_genes = probs.chainSize() - 1;

            // find max J alignment
            segindex_t j_num = clonotype.nJ();
            seq_len_t len = clonotype.sequence().size();
            for (int j_index = 0; j_index < j_num; ++j_index) {
                if (clonotype.getJstart(j_index) < len) {
                    len = clonotype.getJstart(j_index);
                }
            }
            len = clonotype.sequence().size() - len + 1;

            // add J deletions and J gene nodes
            probs.initNode(J_index_dels, j_num, len + 1, 1);
            probs.initNode(J_index_genes,j_num, 1, 1);
            if (full_build) {
                events.initNode(J_index_dels, j_num, len + 1, 1);
                events.initNode(J_index_genes, j_num, 1, 1);
            }

            // compute J deletions
            seq_len_t j_len = 0;
            segindex_t j_gene = 0;
            seq_len_t j_start = 0;

            for (segindex_t j_index = 0; j_index < j_num; ++j_index) {
                j_gene = clonotype.getJ(j_index);
                j_len = _genes->J()[j_gene].sequence.size();
                j_start = clonotype.getJstart(j_index);

                probs(J_index_genes, j_index, 0, 0) = _param_vec->prob_J_gene(j_gene); // probability of choosing this J gene segment
                for (seq_len_t i = 0; i < len + 1; ++i) {
                    if (j_len - len + i >= 0 && j_len - len + i <= clonotype.sequence().size() - j_start) {
                        probs(J_index_dels, j_index, i, 0) = _param_vec->prob_J_del(j_gene, j_len - len + i); // probability of deletions
                    } else {
                        probs(J_index_dels, j_index, i, 0) = 0; // if exceeds length of J gene segment
                    }

                }

                if (full_build) {
                    events(J_index_genes, j_index, 0, 0) = _param_vec->index_J_gene(j_gene);
                    for (seq_len_t i = 0; i < len + 1; ++i) {
                        if (j_len - len - i >= 0) {
                            events(J_index_dels, j_index, i, 0) = _param_vec->index_J_del(j_gene, j_len - len - i);
                        } else {
                            events(J_index_dels, j_index, i, 0) = 0;
                        }
                    }

                    for (seq_len_t i = clonotype.sequence().size() - len + 1; i < clonotype.sequence().size() + 1; ++i) {
                        seq_poses.push_back(i);
                    }
                }
            }
        }


        void buildDiversity(const Clonotype &clonotype,
                ProbMMC &probs,
                EventIndMMC &events,
                vector<seq_len_t> &seq_poses,
                bool full_build) const
        {

        }


        void buildVJinsertions(/*const*/ Clonotype &clonotype,
                ProbMMC &probs,
                EventIndMMC &events,
                vector<seq_len_t> &seq_poses,
                bool full_build) const
        {
            MarkovChain mc(_param_vec->get_iterator(_param_vec->index_VJ_ins_nuc()));

            int v_vertices = probs.nodeColumns(1), j_vertices = probs.nodeRows(3);
            int max_size = _param_vec->max_VJ_ins_len();
            int insertion_len;
            bool good_insertion;

            probs.initNode(VarJoi_INSERTIONS_MATRIX_INDEX, 1, v_vertices, j_vertices);

            if (full_build) {
                events.initNode(VarJoi_INSERTIONS_MATRIX_INDEX, 1, v_vertices, j_vertices);
            }

            for (int v_i = 0; v_i < v_vertices; ++v_i) {
                for (int j_i = v_vertices; j_i < j_vertices + v_vertices; ++j_i) {
                    insertion_len = seq_poses[j_i] - seq_poses[v_i] - 1;
                    good_insertion = (insertion_len >= 0) && (insertion_len < max_size);
                    if (good_insertion) {
                        probs(VarJoi_INSERTIONS_MATRIX_INDEX, 0, v_i, j_i - v_vertices) = mc.nucProbability(clonotype.seq_iterator(seq_poses[v_i]), insertion_len) * _param_vec->prob_VJ_ins_len(insertion_len);
                        if (full_build) {
                            events(VarJoi_INSERTIONS_MATRIX_INDEX, 0, v_i, j_i - v_vertices) = _param_vec->index_VJ_ins_len(insertion_len);
                        }
                    } else {
                        probs(VarJoi_INSERTIONS_MATRIX_INDEX, 0, v_i, j_i - v_vertices) = 0;
                        if (full_build) {
                            events(VarJoi_INSERTIONS_MATRIX_INDEX, 0, v_i, j_i - v_vertices) = 0;
                        }
                    }
                }
            }
        }


        void buildVDinsertions(const Clonotype &clonotype,
                ProbMMC &probs,
                EventIndMMC &events,
                vector<seq_len_t> &seq_poses,
                bool full_build) const
        {
            MarkovChain mc(_param_vec->get_iterator(_param_vec->index_VD_ins_nuc()));
        }


        void buildDJinsertions(const Clonotype &clonotype,
                ProbMMC &probs,
                EventIndMMC &events,
                vector<seq_len_t> &seq_poses,
                bool full_build) const
        {
            MarkovChain mc(_param_vec->get_iterator(_param_vec->index_DJ_ins_nuc()));
        }


        void buildInsertions(Clonotype& clonotype,
                ProbMMC &probs,
                EventIndMMC &events,
                vector<seq_len_t> &seq_poses,
                bool full_build,
                const MarkovChain& mc) {

        }

    };
}

#endif //_MAAGBUILDER_H_
