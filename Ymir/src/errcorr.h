//
// Created by Vadim N. on 27/10/2015.
//

#ifndef YMIR_ERRCORR_H
#define YMIR_ERRCORR_H


#include "clonotype.h"
#include "types.h"
#include "genesegment.h"


namespace ymir {


    class ErrorCorrector {
    public:


        ErrorCorrector(const VDJRecombinationGenes &genes, const SmithWatermanAligner &swa, const SmithWatermanNoGapAligner &swnga)
                : _genes(genes), _swa(swa), _swnga(swnga)
        {

        }


        Clonotype correctAndBuild(const std::string &sequence) const {
            // Smith-Waterman on all alleles of V and J
            SmithWatermanAligner swa();
            // Correct indel errors basing on the max alignment
            // Re-align with Smith-Waterman without indels
            SmithWatermanNoGapAligner swnga();

            // Remove alignments with less than threshold of the max alignment
            // Build clonotype
        }


    protected:

        VDJRecombinationGenes _genes;
        float _threshold;
        AlignmentEventScore _v_score, _d_score, _j_score;
        SmithWatermanAligner _swa;
        SmithWatermanNoGapAligner _swnga;


        ErrorCorrector() {}


        void smith_waterman(const std::string &sequence) const {

        }

    };
}

#endif //YMIR_ERRCORR_H
