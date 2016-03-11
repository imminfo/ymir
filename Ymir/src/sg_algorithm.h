//
// Created by Vadim N. on 11/03/2016.
//

#ifndef YMIR_SG_ALGORITHM_H
#define YMIR_SG_ALGORITHM_H


#include "em_algorithm.h"


namespace ymir {

    class BlockEMAlgorithm;

    /**
    * \class BlockEMAlgorithm
    *
    * \brief Implementation of the online EM-algorithm for statistical inference of assembling model parameters.
    */
    class BlockEMAlgorithm : public EMAlgorithm {
    public:

        virtual bool statisticalInference(const ClonesetView& repertoire,
                                          ProbabilisticAssemblingModel & model,
                                          const AlgorithmParameters& algo_param = AlgorithmParameters().set("niter", 10).set("block.size", 2000),
                                          ErrorMode error_mode = NO_ERRORS) const
        {
            // shuffle input data at each step
            // subvec

            std::vector<size_t> indices;
            size_t start_i;
            size_t block_size = algo_param["block.size"].asUInt();
//            prob_t step = algo_param["block.size"].asDouble();
            ModelParameterVector new_param_vec = model.event_probabilities();
            std::vector<bool> changed(new_param_vec.size(), false);

            ClonesetView rep_nonc = repertoire.noncoding().shuffle();
            auto maag_rep = model.buildGraphs(rep_nonc, SAVE_METADATA, error_mode, NUCLEOTIDE, true);

            std::vector<prob_t> prob_vec(maag_rep.size(), 0);
            prob_t prev_ll = 0;

            MAAGForwardBackwardAlgorithm fb;
            for (size_t iter = 1; iter <= algo_param["niter"].asUInt(); ++iter) {

                new_param_vec.fill(0);

                if (start_i + block_size > maag_rep.size()) {
                    start_i = 0;
                    this->getRandomIndices(indices, maag_rep.size());
                }

                for (size_t maag_i = start_i; maag_i < std::min(maag_rep.size(), start_i + block_size); ++maag_i) {
                    // compute marginal probabilities for this block
                    // and update the temporary model parameter vector
                    this->updateTempVec(fb, maag_rep[indices[maag_i]], new_param_vec, changed, error_mode);
                }

                this->updateModel(model, new_param_vec, maag_rep, prob_vec, prev_ll, changed, error_mode);

                // update starting index
                start_i += block_size;
            }
        }

    protected:

        vector<size_t> getRandomIndices(vector<size_t> &indices, size_t size) const {

        }


        bool updateTempVec(MAAGForwardBackwardAlgorithm &fb,
                           MAAG &maag,
                           ModelParameterVector &new_param_vec,
                           vector<bool> &changed,
                           ErrorMode error_mode) const
        {
            if (!fb.process(maag, error_mode)) {
                return false;
            }

            while (!fb.is_empty()) {
                event_pair_t ep = fb.nextEvent();
                new_param_vec[ep.first] += ep.second;
                changed[ep.first] = true;
            }

            if (error_mode) { new_param_vec.set_error_prob(new_param_vec.error_prob() + fb.err_prob()); }

            if (maag.is_vj()) {
                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)] += fb.VJ_nuc_probs()[0];
                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 1)] += fb.VJ_nuc_probs()[1];
                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 2)] += fb.VJ_nuc_probs()[2];
                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 3)] += fb.VJ_nuc_probs()[3];

                changed[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)] = true;
                changed[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 1)] = true;
                changed[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 2)] = true;
                changed[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 3)] = true;
            } else {
                int k = 0;
                for (auto prev_nuc = 0; prev_nuc < 4; ++prev_nuc) {
                    for (auto next_nuc = 0; next_nuc < 4; ++next_nuc, ++k) {
                        new_param_vec[new_param_vec.event_index(VDJ_VAR_DIV_INS_NUC, prev_nuc, next_nuc)] += fb.VD_nuc_probs()[k];
                        changed[new_param_vec.event_index(VDJ_VAR_DIV_INS_NUC, prev_nuc, next_nuc)] = true;
                    }
                }

                k = 0;
                for (auto prev_nuc = 0; prev_nuc < 4; ++prev_nuc) {
                    for (auto next_nuc = 0; next_nuc < 4; ++next_nuc, ++k) {
                        new_param_vec[new_param_vec.event_index(VDJ_DIV_JOI_INS_NUC, prev_nuc, next_nuc)] += fb.DJ_nuc_probs()[k];
                        changed[new_param_vec.event_index(VDJ_DIV_JOI_INS_NUC, prev_nuc, next_nuc)] = true;
                    }
                }
            }
        }


        void updateModel(ProbabilisticAssemblingModel &model,
                         ModelParameterVector &new_param_vec,
                         MAAGRepertoire &maag_rep,
                         vector<prob_t> &prob_vec,
                         prob_t &prev_ll,
                         vector<bool> &changed,
                         ErrorMode error_mode) const
        {
            if (error_mode) { new_param_vec.set_error_prob(new_param_vec.error_prob() / (maag_rep.size())); }

            new_param_vec.normaliseEventFamilies();

            for (size_t i = 0; i < new_param_vec.size(); ++i) {
                if (!changed[i]) {
                    new_param_vec[i] = model.event_probabilities()[i];
                }
            }

            model.updateModelParameterVector(new_param_vec);
            model.updateEventProbabilities(&maag_rep);

            for (size_t i = 0; i < maag_rep.size(); ++i) {
                prob_vec[i] = maag_rep[i].fullProbability();
            }
            prob_summary(prob_vec, prev_ll);
            prev_ll = loglikelihood(prob_vec);

            changed.clear();
            changed.resize(new_param_vec.size(), false);
        }


    };

}

#endif //YMIR_SG_ALGORITHM_H
