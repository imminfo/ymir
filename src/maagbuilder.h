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


#include "maag.h"


namespace ymir {

    class MAAGBuilder;


    /**
    * \class MAAGBuilder
    */
    class MAAGBuilder : protected MAAG {

    public:


        MAAGBuilder(const ModelParameterVector &param_vec, const VDJRecombinationGenes &genes) {
            _param_vec = new ModelParameterVector(param_vec);
            _genes = new VDJRecombinationGenes(genes);
        }


        virtual ~MAAGBuilder() {
            if (_param_vec) { delete _param_vec; }
            if (_genes) { delete _genes; }
        }


        MAAG build(const Clonotype &clonotype, bool full_build = false) const {
            ProbMMC probs;
            EventIndMMC events;
            vector<seq_len_t> seq_poses;
            seq_poses.reserve(DEFAULT_SEQ_POSES_RESERVE);


            if (!clonotype.is_vdj()) {
                probs.resize(5);
                if (full_build) { events.resize(5); }
            } else {
                probs.resize(7);
                if (full_build) { events.resize(7); }
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
            for (int i = 0; i < seq_poses.size(); ++i) { seq_poses_arr[i] = seq_poses[i]; }

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


        MAAGBuilder() {}


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
            probs.initNode(0, v_num, 1, 1);
            probs.initNode(1, v_num, 1, len + 1);
            if (full_build) {
                events.initNode(0, v_num, 1, 1);
                events.initNode(1, v_num, 1, len + 1);
            }
            for (segindex_t v_index = 0; v_index < v_num; ++v_index) {
                v_gene = clonotype.getV(v_index);
                v_len = _genes->V()[v_gene].sequence.size();

                probs(0, v_index, 0, 0) = _param_vec->prob_V_gene(v_gene); // probability of choosing this V gene segment
                for (seq_len_t i = 0; i < len + 1; ++i) {
                    if (v_len - i >= 0) {
                        probs(1, v_index, 0, i) = _param_vec->prob_V_del(v_gene, v_len - i); // probability of deletions
                    } else {
                        probs(1, v_index, 0, i) = 0; // if exceeds length of V gene segment
                    }

                }

                if (full_build) {
                    events(0, v_index, 0, 0) = _param_vec->index_V_gene(v_gene);
                    for (seq_len_t i = 0; i < len + 1; ++i) {
                        if (v_len - i >= 0) {
                            events(1, v_index, 0, i) = _param_vec->index_V_del(v_gene, v_len - i);
                        } else {
                            events(1, v_index, 0, i) = _param_vec->index_V_del(v_gene, 0);
                        }
                    }

                    for (seq_len_t i = 0; i < len + 1; ++i) {
                        seq_poses.push_back(i);
                    }
                }
            }

        }


        void buildJoining(const Clonotype &clonotype,
                ProbMMC &probs,
                EventIndMMC &events,
                vector<seq_len_t> &seq_poses,
                bool full_build) const
        {
            int J_index_dels, J_index_genes;
            if (!clonotype.is_vdj()) {
                J_index_dels = 3;
                J_index_genes = 4;
            } else {
                J_index_dels = 5;
                J_index_genes = 6;
            }

            // find max J alignment
            // WRONG BECAUSE WE HAVE J START NOT J ALIGNEMNT LENGTH
            segindex_t j_num = clonotype.nJ();
            seq_len_t len = 0;
            for (int j_index = 0; j_index < j_num; ++j_index) {
                if (clonotype.getVend(j_index) > len) {
                    len = clonotype.getJstart(j_index);
                }
            }

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

            for (segindex_t j_index = 0; j_index < j_num; ++j_index) {
                j_gene = clonotype.getJ(j_index);
                j_len = _genes->J()[j_gene].sequence.size();

                int J_DELS = 3;
                int J_INDEX = 4;

                probs(J_INDEX, j_index, 0, 0) = _param_vec->prob_J_gene(j_gene); // probability of choosing this J gene segment
                for (seq_len_t i = len; i >= 0; --i) {
                    if (j_len - i >= 0) {
                        probs(J_DELS, j_index, i, 0) = _param_vec->prob_J_del(j_gene, j_len - i); // probability of deletions
                    } else {
                        probs(J_DELS, j_index, i, 0) = 0; // if exceeds length of J gene segment
                    }

                }

                if (full_build) {
                    events(J_INDEX, j_index, 0, 0) = _param_vec->index_J_gene(j_gene);
                    for (seq_len_t i = len; i >= 0; --i) {
                        if (j_len - i >= 0) {
                            events(J_DELS, j_index, i, 0) = _param_vec->index_J_del(j_gene, j_len - i);
                        } else {
                            events(J_DELS, j_index, i, 0) = 0;
                        }
                    }

                    for (seq_len_t i = 0; i < len + 1; ++i) {
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


        void buildVJinsertions(const Clonotype &clonotype,
                ProbMMC &probs,
                EventIndMMC &events,
                vector<seq_len_t> &seq_poses,
                bool full_build) const
        {
            MarkovChain mc(_param_vec->get_iterator(_param_vec->index_VJ_ins_nuc()));
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

    };
}

#endif //_MAAGBUILDER_H_
