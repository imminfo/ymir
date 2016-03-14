//
// Created by Vadim N. on 11/03/2016.
//

#ifndef YMIR_SG_ALGORITHM_H
#define YMIR_SG_ALGORITHM_H


#include "em_algorithm.h"


namespace ymir {

    class SGAlgorithm;


    /**
    * \class BlockEMAlgorithm
    *
    * \brief Implementation of the online EM-algorithm for statistical inference of assembling model parameters.
    */
    class SGAlgorithm : public EMAlgorithm {
    public:

        virtual std::vector<prob_t> statisticalInference(const ClonesetView& repertoire,
                                                         ProbabilisticAssemblingModel & model,
                                                         const AlgorithmParameters& algo_param = AlgorithmParameters()
                                                                 .set("niter", 10)
                                                                 .set("block.size", 5000)
                                                                 .set("alpha", .6)
                                                                 .set("prebuild", false)
                                                                 .set("sample", 50000),
                                                         ErrorMode error_mode = NO_ERRORS) const
        {
            // shuffle input data at each step
            // subvec -4008648.1

            cout << "Statistical inference on a PAM:\t" << model.name() << endl;
            cout << "\tOnline EM-algorithm.";
            if (error_mode == COMPUTE_ERRORS) {
                cout << "\t(with sequence errors)";
            }
            std::cout << std::endl;


            if (!algo_param.check("niter") && !algo_param.check("block.size") && !algo_param.check("alpha") && !algo_param.check("sample")) {
                return std::vector<prob_t>();
            }

            std::cout << "#iterations:\t" << (size_t) algo_param["niter"].asUInt() << std::endl;
            std::cout << "block size:\t" << (size_t) algo_param["block.size"].asUInt() << std::endl;
            std::cout << "alpha:\t" << (size_t) algo_param["block.size"].asUInt() << std::endl;
            std::cout << "prebuild:\t" << (size_t) algo_param["prebuild"].asBool() << std::endl;

            std::vector<prob_t> logLvec;

            size_t sample = algo_param["sample"].asUInt();
            ClonesetView rep_nonc = repertoire.noncoding().sample(sample);
            cout << "Number of noncoding clonotypes:\t" << (size_t) rep_nonc.size() << endl;


            size_t start_i = rep_nonc.size();
            size_t block_size = algo_param["block.size"].asUInt();
            prob_t alpha = algo_param["alpha"].asDouble(); // step(k) = (k + 2)^(-alpha), .5 < alpha <= 1
            ModelParameterVector new_param_vec = model.event_probabilities();
            new_param_vec.fill(1);
//            new_param_vec.set_error_prob(.003);
            new_param_vec.normaliseEventFamilies();
            model.updateModelParameterVector(new_param_vec);

            std::vector<bool> changed(new_param_vec.size(), false);

            auto maag_rep = model.buildGraphs(rep_nonc, SAVE_METADATA, error_mode, NUCLEOTIDE, true);

            std::vector<prob_t> prob_vec(maag_rep.size(), 0);
            prob_t prev_ll = 0;

            cout << "Computing full assembling probabilities..." << endl;
            vector<bool> good_clonotypes;
            size_t removed, zero_prob, no_alignments;
            this->filterOut(rep_nonc, maag_rep, prob_vec, good_clonotypes, removed, zero_prob, no_alignments);

            std::vector<size_t> indices;
            for (size_t i = 0; i < rep_nonc.size(); ++i) {
                if (good_clonotypes[i]) {
                    indices.push_back(i);
                }
            }
            std::random_shuffle(indices.begin(), indices.end());

            cout << endl << "Initial data summary:" << endl;
            prob_summary(prob_vec);
            std::cout << model.event_probabilities().error_prob() << std::endl;
            prev_ll = loglikelihood(prob_vec);
            logLvec.push_back(prev_ll);

            MAAGForwardBackwardAlgorithm fb;
            for (size_t iter = 1; iter <= algo_param["niter"].asUInt(); ++iter) {
                if (start_i + block_size > indices.size()) {
                    start_i = 0;
                    std::random_shuffle(indices.begin(), indices.end());
                } else {
                    start_i += block_size;
                }

                std::cout << "=======================" << std::endl
                << "Iteration: " << (size_t) iter << " Block: [" << (int) start_i << ":" << (int) std::min(indices.size() - 1, start_i + block_size - 1)
                << "]" << std::endl << "=======================" << std::endl;

                new_param_vec.fill(0);

                for (size_t maag_i = start_i; maag_i < std::min(indices.size(), start_i + block_size); ++maag_i) {
                    // compute marginal probabilities for this block
                    // and update the temporary model parameter vector
                    this->updateTempVec(fb, maag_rep[indices[maag_i]], new_param_vec, changed, error_mode);
                }

                this->updateModel(model, new_param_vec, maag_rep, prob_vec, prev_ll, changed, exp(-alpha * log(iter + 2)), error_mode);

                logLvec.push_back(prev_ll);
            }

            return logLvec;
        }

    protected:

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
//                if (std::isnan(fb.VJ_nuc_probs()[0])) { cout << "A" << endl; check_and_throw(std::isnan(fb.VJ_nuc_probs()[0]), "123"); }
//                if (std::isnan(fb.VJ_nuc_probs()[1])) {cout << "C" << endl; check_and_throw(std::isnan(fb.VJ_nuc_probs()[1]), "123"); }
//                if (std::isnan(fb.VJ_nuc_probs()[2])) {cout << "G" << endl; check_and_throw(std::isnan(fb.VJ_nuc_probs()[2]), "123"); }
//                if (std::isnan(fb.VJ_nuc_probs()[3])) {cout << "T" << endl; check_and_throw(std::isnan(fb.VJ_nuc_probs()[3]), "123"); }
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
                         prob_t step_k,
                         ErrorMode error_mode) const
        {
            if (error_mode) { new_param_vec.set_error_prob(step_k * model.event_probabilities().error_prob() + (1 - step_k) * new_param_vec.error_prob()); }

            new_param_vec.normaliseEventFamilies();

            for (size_t i = 0; i < new_param_vec.size(); ++i) {
                if (!changed[i]) {
                    new_param_vec[i] = model.event_probabilities()[i];
                } else {
                    new_param_vec[i] = step_k * model.event_probabilities()[i] + (1 - step_k) * new_param_vec[i];
                }
            }

            new_param_vec.normaliseEventFamilies();

            model.updateModelParameterVector(new_param_vec);
            model.updateEventProbabilities(&maag_rep, false);

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
