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


        MAAG build(const Clonotype &clonotype, bool full_build = false) {
            ProbMMC probs;
            EventIndMMC events;
            vector<seq_len_t> seq_poses;
            seq_poses.reserve(300);

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
            probs.addNode(v_num, 1, 1);
            probs.addNode(v_num, 1, len + 1);
            if (full_build) {
                events.addNode(v_num, 1, 1);
                events.addNode(v_num, 1, len + 1);
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


            // find max J alignment
            // WRONG BECAUSE WE HAVE J START NOT J ALIGNEMNT LENGTH
            segindex_t j_num = clonotype.nJ();
            len = 0;
            for (int j_index = 0; j_index < j_num; ++j_index) {
                if (clonotype.getVend(j_index) > len) {
                    len = clonotype.getJstart(j_index);
                }
            }


            // VJ - recombination
            if (!clonotype.is_vdj()) {

                // add VJ insertions node
                probs.addNode(1, probs.nodeColumns(1), len + 1);
                if (full_build) {
                    events.addNode(1, probs.nodeColumns(1), len + 1);
                }

                // add J deletions and J gene nodes
                probs.addNode(j_num, len + 1, 1);
                probs.addNode(j_num, 1, 1);
                if (full_build) {
                    events.addNode(j_num, len + 1, 1);
                    events.addNode(j_num, 1, 1);
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



                // compute VJ insertions
                MarkovChain mc(_param_vec->get_iterator(_param_vec->index_VJ_ins_nuc()));


            }
            // VDJ - recombination
            else {
                // get all possible D alignments


                // compute VD insertions
                MarkovChain mc_vd(_param_vec->get_iterator(_param_vec->index_VD_ins_nuc()));


                // compute J deletions


                // compute DJ insertions
                MarkovChain mc_dj(_param_vec->get_iterator(_param_vec->index_DJ_ins_nuc()));

            }

            if (full_build) {
                return MAAG(probs, events, clonotype.sequence(), seq_poses, seq_poses.size());
            } else {
                return MAAG(probs);
            }
        }


//        MAAGRepertoire build(const Cloneset &cloneset, bool full_build = false) {
//            // ???
//        }

    // vector<numeric> buildAndCompute(const Cloneset &cloneset)


    protected:

        ModelParameterVector *_param_vec;  // or just copy it?
        VDJRecombinationGenes *_genes; // copy this too?


        MAAGBuilder() {}

    };
}

#endif //_MAAGBUILDER_H_
