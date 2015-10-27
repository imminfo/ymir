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


        struct Score {
            float match, mismatch, ins, del;

            Score()
                    : match(1), mismatch(-1), ins(-1), del(-1)
            {}

            Score(float match_, float mismatch_, float ins_, float del_)
                    : match(match_), mismatch(mismatch_), ins(ins_), del(del_)
            {}
        };


        ErrorCorrector(const VDJRecombinationGenes &genes,
                       const Score &v_score = Score(1, -1, -1, -1),
                       const Score &d_score = Score(1, -1, -1, -1),
                       const Score &j_score = Score(1, -1, -1, -1))
                : _genes(genes), _v_score(v_score), _d_score(d_score), _j_score(j_score)
        {

        }


        Clonotype correctAndBuild(const std::string &sequence) const {
            // Smith-Waterman on all alleles of V and J
            // Correct indel errors basing on the max alignment
            // Re-align with Smith-Waterman without indels
            // Remove alignments with less than threshold of the max alignment
            // Build clonotype
        }


    protected:

        VDJRecombinationGenes _genes;
        float _threshold;
        Score _v_score, _d_score, _j_score;


        ErrorCorrector() {}


        void smith_waterman(const std::string &sequence) const {

        }

    };
}

#endif //YMIR_ERRCORR_H
