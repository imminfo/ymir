/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vadim dot nazarov at mailbox dot com>
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


    class MAAGBuilder : protected MAAG {

    protected:

        /**
        *
        */
        typedef tuple<segindex_t, prob_t, eventind_t, segindex_t> gene_info;

    public:


        MAAGBuilder() : MAAG() {

        }


        virtual ~MAAGBuilder() {
        }


        MAAG build() {

        }


        void reset() { delete _maag; }


        void addVariableGene(segindex_t v_index
                , eventind_t v_event
                , prob_t v_prob) {

        }


        void addVariableDeletion(eventind_t v_order
                , eventind_t event_index
                , prob_t event_prob) {

        }

        void addJoiningGene() {}

        void addDiversityGene() {}


    protected:

        MAAG *_maag;
        vector<gene_info> _vgenes, _jgenes, _dgenes;
    };
}

#endif //_MAAGBUILDER_H_
