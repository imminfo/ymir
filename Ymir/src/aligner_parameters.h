//
// Created by Vadim N. on 11/03/2016.
//

#ifndef YMIR_ALIGNER_PARAMETERS_H
#define YMIR_ALIGNER_PARAMETERS_H


#include "types.h"


namespace ymir {


    struct AlignmentEventScore;
    struct VDJAlignmentEventScore;
    struct VDJAlignmentScoreThreshold;
    struct VDJAlignerParameters;


    /**
     *
     */
    struct AlignmentEventScore {

        /**
         *
         */
        AlignmentEventScore(alignment_score_t match_,
                            alignment_score_t mism_,
                            alignment_score_t indel_)
            : match(match_),
              mism(mism_),
              indel(indel_)
        {
        }


        alignment_score_t match, mism, indel;

    };


    struct VDJAlignmentEventScore {

        static const alignment_score_t default_v_match = 2;
        static const alignment_score_t default_v_mism = -3;
        static const alignment_score_t default_v_indel = -3;

        static const alignment_score_t default_d_match = 2;
        static const alignment_score_t default_d_mism = -3;
        static const alignment_score_t default_d_indel = -3;

        static const alignment_score_t default_j_match = 2;
        static const alignment_score_t default_j_mism = -3;
        static const alignment_score_t default_j_indel = -3;


        VDJAlignmentEventScore()
            : v_score(default_v_match, default_v_mism, default_v_indel),
              d_score(default_d_match, default_d_mism, default_d_indel),
              j_score(default_j_match, default_j_mism, default_j_indel)
        {
        }


        VDJAlignmentEventScore(const AlignmentEventScore &v,
                               const AlignmentEventScore &j)
            : v_score(v),
              d_score(0, 0, 0),
              j_score(j)
        {
        }

        VDJAlignmentEventScore(const AlignmentEventScore &v,
                               const AlignmentEventScore &d,
                               const AlignmentEventScore &j)
            : v_score(v),
              d_score(d),
              j_score(j)
        {
        }

        AlignmentEventScore v_score, d_score, j_score;

    };


    struct VDJAlignmentScoreThreshold {

        static const alignment_score_t default_v_threshold = 20;
        static const alignment_score_t default_d_threshold = 20;
        static const alignment_score_t default_j_threshold = 20;

        VDJAlignmentScoreThreshold()
            : v_threshold(default_v_threshold),
              d_threshold(default_d_threshold),
              j_threshold(default_j_threshold)
        {
        }


        VDJAlignmentScoreThreshold(alignment_score_t v, alignment_score_t j)
            : v_threshold(v),
              d_threshold(0),
              j_threshold(j)
        {
        }


        VDJAlignmentScoreThreshold(alignment_score_t v, alignment_score_t d, alignment_score_t j)
            : v_threshold(v),
              d_threshold(d),
              j_threshold(j)
        {
        }


        alignment_score_t v_threshold, d_threshold, j_threshold;

    };


    /**
     * \struct VDJAlignerParameters
     */
    struct VDJAlignerParameters {

        static const seq_len_t default_minlen = 3;


        VDJAlignerParameters()
            : min_D_len(default_minlen) {
        }


        VDJAlignerParameters(seq_len_t minlen)
            : min_D_len(minlen) {
        }



        VDJAlignerParameters(seq_len_t minlen,
                             VDJAlignmentEventScore score_,
                             VDJAlignmentScoreThreshold thr_)
            : min_D_len(minlen),
              threshold(thr_),
              score(score_)
        {
        }


        seq_len_t min_D_len; ///
        VDJAlignmentEventScore score;
        VDJAlignmentScoreThreshold threshold;

    };

}

#endif //YMIR_ALIGNER_PARAMETERS_H
