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


#include <AVFoundation/AVFoundation.h>
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

            // find max alignment
            seq_len_t len = 0;
            for (int v_index = 0; v_index < clonotype.nV(); ++v_index) {
                if (clonotype.getVend(v_index) > len) {
                    len = clonotype.getVend(v_index);
                }
            }

            // compute V deletions
            segindex_t v_num = clonotype.nV();
            seq_len_t v_len = 0;
            segindex_t v_gene = 0;
            probs.addNode(v_num, 1, 1);
            probs.addNode(v_num, 1, len + 1);
            if (full_build) {
                events.addNode(v_num, 1, 1);
                events.addNode(v_num, 1, len + 1);
            }
            for (segindex_t v_index = 0; v_index < clonotype.nV(); ++v_index) {
                v_gene = clonotype.getV(v_index);
                v_len = _genes->V()[v_gene].sequence.size();

                probs(0, v_index, 0, 0) = _param_vec->prob_V_gene(v_gene); // probability of choosing this V gene segment
                for (seq_len_t i = 0; i < len + 1; ++i) {
                    probs(1, v_index, 0, i) = _param_vec->prob_V_del(v_gene, v_len - i); // probability of deletions
                }

                if (full_build) {
                    events(0, v_index, 0, 0) = _param_vec->index_V_gene(v_gene);
                    for (seq_len_t i = 0; i < len + 1; ++i) {
                        events(1, v_index, 0, i) = _param_vec->index_V_del(v_gene, v_len - i);
                    }

                    for (seq_len_t i = 0; i < len + 1; ++i) {
                        seq_poses.push_back(i);
                    }
                }
            }

            // compute probability of all insertions
            // VJ - recombination
            if (!clonotype.is_vdj()) {
                // compute J deletions


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
                return MAAG(probs, events, clonotype.sequence(), seq_poses, n_poses);
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


        MAAGBuilder() : MAAG() {}

    };
}

#endif //_MAAGBUILDER_H_
