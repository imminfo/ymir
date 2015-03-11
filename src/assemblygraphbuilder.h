/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vdm dot nazarov at gmail dot com>
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

#ifndef _ASSEMBLY_GRAPH_BUILDER_H
#define _ASSEMBLY_GRAPH_BUILDER_H


#include "assemblygraph.h"
#include "modelparametervector.h"


namespace ymir {

    /*
    Builder:

    nucleotide sequence builder

    amino acid sequence builder

    build()
        Clone, ModelParameterVector
            -> AssemblyGraph
        ClonalRepertoireView, ModelParameterVector
            -> AssemblyGraphRepertoire

    replaceEventProbabilities()
        &AssemblyGraph, ModelParameterVector
            -> bool
        &AssemblyGraphRepertoire, ModelParameterVector
            -> bool
    */


    class AssemblyGraphBuilder : protected AssemblyGraph {

    public:


        AssemblyGraphBuilder() : AssemblyGraph() {

        }


        virtual ~AssemblyGraphBuilder() {
        }


        static AssemblyGraph buildGraph(const CloneMetadata& clone, const ModelParameterVector& param_vec) const {

        }


        static vector<AssemblyGraph> buildGraph() const {

        }


        // build MultiAlignmentAssemblyGraph
        // buildMAGraph


        // build AminoAcidMotifAssemblyGraph
        // buildAAMotifGraph

    protected:

    };
}

#endif