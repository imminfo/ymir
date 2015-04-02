//
// Created by Vadim N. on 24/03/2015.
//

#ifndef _YMIR_CLONOTYPEASSEMBLER_H_
#define _YMIR_CLONOTYPEASSEMBLER_H_


#include "clonotype.h"


namespace ymir {

    class ClonotypeAssembler {
    public:

        /**
         *
         */
        ClonotypeAssembler( /* ModelParameterVector, VDJRecombinationGenes */ ) { }

    protected:


    };

//    /**
//        * \class SequenceGenerator
//        *
//        * \brief Generator of artificial sequences.
//        */
//    struct SequenceGenerator { // markov chain
//
//    protected:
//
//        struct GeneSegmentNode {
//            // shared ptrs to next nodes w/ probability of choosing
//        };
//
//
//        struct DeletionNode {
//            // shared ptrs to next nodes with insertion probabilities
//        };
//
//
//        struct InsertionNode {
//
//        };
//
//    public:
//
//        SequenceGenerator() {
//
//        }
//
//
//        SequenceGenerator& addGeneMatrix(const event_matrix_t& mat) {
//
//            return *this;
//        }
//
//
//        SequenceGenerator& addDelMatrix(const event_matrix_t& mat) {
//
//            return *this;
//        }
//
//
//        /**
//        * \brief Add common insertions probability matrix for all previous gene segments.
//        */
//        SequenceGenerator& addInsVector(const vector<prob_t> &vec, const MarkovChain& chain) {
//
//            return *this;
//        }
//
//
////            ClonalRepertoire generate(size_t count = 1) const {
//        // to generate a Clone just make a path through a graph with events
////            }
//
//    protected:
//
//    };
}

#endif //_YMIR_CLONOTYPEASSEMBLER_H_
