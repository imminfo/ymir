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

#ifndef _CLONOTYPE_BUILDER_H
#define _CLONOTYPE_BUILDER_H


#include "clonotype.h"
#include "vdj_alignment_builder.h"


namespace ymir {

    #define CLONOTYPEBUILDER_VSEG_DEFAULT_RESERVE_SIZE 9
    #define CLONOTYPEBUILDER_JSEG_DEFAULT_RESERVE_SIZE 6
    #define CLONOTYPEBUILDER_DSEG_DEFAULT_RESERVE_SIZE 30


    /**
    * \class ClonotypeBuilder
    */
    class ClonotypeBuilder : protected Clonotype, public VDJAlignmentBuilder {
    public:

        ClonotypeBuilder() : VDJAlignmentBuilder()
        {
        }


        virtual ~ClonotypeBuilder() {

        }


        /**
         * \brief Build clone alignment structure with stored information.
         *
         * \return Pointer to the newly created ClonotypeAlignment object.
        */
        Clonotype buildClonotype() {
//            ClonotypePtr cla(new Clonotype(_sequence, _seq_type, _recomb, this->buildAlignment()));
//            return std::move(cla);
            return Clonotype(_sequence, _seq_type, _recomb, this->buildAlignment());
        }


        /**
         *
         */
        ///@{
        ClonotypeBuilder& setSequence(const std::string& seq) { this->_sequence = seq; return *this; }

        ClonotypeBuilder& setSequenceType(SequenceType seq_type) { _seq_type = seq_type; return *this; }

        ClonotypeBuilder& setNucleotideSeq() { _seq_type = NUCLEOTIDE; return *this; }

        ClonotypeBuilder& setAminoAcidSeq() { _seq_type = AMINOACID; return *this; }

        ClonotypeBuilder& setRecombination(Recombination recomb) { _recomb = recomb; return *this; }
        ///@}

    protected:

    };

}

#endif