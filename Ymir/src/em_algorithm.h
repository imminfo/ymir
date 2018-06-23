//
// Created by Vadim N. on 11/03/2016.
//

#ifndef YMIR_EM_ALGORITHM_H
#define YMIR_EM_ALGORITHM_H


#include "statisticalinferencealgorithm.h"


namespace ymir {

    class EMAlgorithm;


    /**
     * \class EMAlgorithm
     *
     * \brief Implementation of the EM-algorithm for statistical inference of assembling model parameters.
     * Classic version described in (Murugan et al 2012)
     */
    class EMAlgorithm : public StatisticalInferenceAlgorithm {
    public:

        /**
         *
         */
        virtual std::vector<prob_t> statisticalInference(const ClonesetViewNuc &repertoire,
                                                         ProbabilisticAssemblingModel &model,
                                                         const AlgorithmParameters &algo_param = AlgorithmParameters().set("niter", 10).set("sample", 50000).set("memory-safe", false),
                                                         ErrorMode error_mode = NO_ERRORS) const {
            cout << "Statistical inference on a PAM:\t" << model.name() << endl;
            cout << "\nBatch EM-algorithm.";
            if (error_mode == COMPUTE_ERRORS) {
                cout << "\t(with sequence errors)";
            }
            std::cout << std::endl;

            if (!algo_param.check("niter") && !algo_param.check("sample") && !algo_param.check("memory-safe")) {
                return std::vector<prob_t>();
            }

            cout << "\t#iterations: "
                 << (int) algo_param["niter"]
                 << std::endl;

            size_t sample = algo_param["sample"];
            if (sample) {
                std::cout << "\tsample size (doesn't work yet): "
                          << (int) sample
                          << std::endl;
            }

            bool memory_safe = algo_param["memory-safe"];
            if (memory_safe) {
                std::cout << "\tworking in the memory-safe mode: rebuild graphs at each step"
                          << std::endl;
            } else {
                std::cout << "\tworking in the speed-wise mode: build graphs only once"
                          << std::endl;
            }

            ClonesetViewNuc rep_nonc = repertoire.noncoding();
            cout << "Number of noncoding clonotypes:\t" << (size_t) rep_nonc.size() << endl;

            std::vector<prob_t> logLvec;


            ModelParameterVector new_param_vec = model.event_probabilities();
            new_param_vec.fill(1);
            // Error probability
            new_param_vec.set_error_prob(.0003);

            new_param_vec.normaliseEventFamilies();
            model.updateModelParameterVector(new_param_vec);

            MAAGNucRepertoire maag_rep;
            if (!memory_safe) {
                maag_rep = model.buildGraphs(rep_nonc, SAVE_METADATA, error_mode, true);
            }

            vector<prob_t> prob_vec;

            prob_t prev_ll = 0, cur_ll = 0;

            vector<bool> good_clonotypes;
            size_t removed, zero_prob, no_alignments;
            this->filterOut(rep_nonc, model, maag_rep, prob_vec, good_clonotypes, removed, zero_prob, no_alignments, error_mode, memory_safe);

            cout << endl << "Initial data summary:" << endl;
            prob_summary(prob_vec);
            std::cout << model.event_probabilities().error_prob() << std::endl;
            prev_ll = loglikelihood(prob_vec);
            logLvec.push_back(prev_ll);

            std::chrono::system_clock::time_point tp1, tp2;
            tp1 = std::chrono::system_clock::now();
            for (size_t iter = 1; iter <= algo_param["niter"]; ++iter) {
                if (iter == 1) {
                    cout << endl << "Iteration:\t1 / " << (size_t) algo_param["niter"] << endl;
                } else {
                    cout << endl << "Iteration:\t" << (size_t) iter << " / " << (size_t) algo_param["niter"]
                    << print_time(tp1, algo_param["niter"], iter)
                    << endl;
                }

                new_param_vec.fill(0);


                //
                // Memory safe inference - rebuild MAAGs at each step
                //
                std::cout << "Infer parameters..." << std::endl;
                if (memory_safe) {
#ifdef USE_OMP
                    auto max_thrs = omp_get_max_threads();

                    std::vector<size_t> blocks;
                    size_t block_step = std::min(rep_nonc.size() / max_thrs + 1, rep_nonc.size());
                    blocks.push_back(0);
                    blocks.push_back(block_step);
                    for (auto i = 1; i < max_thrs; ++i) {
                        blocks.push_back(*--blocks.end());
                        blocks.push_back(std::min(*--blocks.end() + block_step, rep_nonc.size()));
                    }

                    #pragma omp parallel
                    {
                        vector<MAAGForwardBackwardAlgorithm> fb(max_thrs);
                        vector<ModelParameterVector> local_param_vec;
                        for (auto i = 0; i < max_thrs; ++i) { local_param_vec.push_back(new_param_vec); }

                        int tid = omp_get_thread_num();

                        size_t start_i = blocks[tid*2],
                               end_i = blocks[tid*2 + 1];

                        for (size_t i = start_i; i < end_i; ++i) {
                            if (good_clonotypes[i]) {
                                this->updateTempVec(fb[tid], model.buildGraphs(rep_nonc[i], SAVE_METADATA, error_mode), local_param_vec[tid], error_mode);
                            }
                        }

                        for (size_t i = 0; i < new_param_vec.size(); ++i) {
                            #pragma omp atomic
                            new_param_vec[i] += local_param_vec[tid][i];
                        }
                    }
#else
                    MAAGForwardBackwardAlgorithm fb;
                    tp1 = std::chrono::system_clock::now();
                    for (size_t i = 0; i < rep_nonc.size(); ++i) {
                        if ((i+1) % 25000 == 0) {
                            cout << "Processed " << (int) i << " / " << (int) rep_nonc.size() << " MAAGs.\t" << endl;
                        }
                        if (good_clonotypes[i]) {
                            if(!this->updateTempVec(fb, model.buildGraphs(rep_nonc[i], SAVE_METADATA, error_mode), new_param_vec, error_mode)) {
                                cout << "bad maag:\t" << (int) i << endl;
                            }
                        }
                    }
#endif
                }

                //
                // Speed-wise inference - build MAAGs only once
                //
                else {
#ifdef USE_OMP
                    auto max_thrs = omp_get_max_threads();

                    std::vector<size_t> blocks;
                    size_t block_step = std::min(rep_nonc.size() / max_thrs + 1, rep_nonc.size());
                    blocks.push_back(0);
                    blocks.push_back(block_step);
                    for (auto i = 1; i < max_thrs; ++i) {
                        blocks.push_back(*--blocks.end());
                        blocks.push_back(std::min(*--blocks.end() + block_step, rep_nonc.size()));
                    }

                    #pragma omp parallel
                    {
                        vector<MAAGForwardBackwardAlgorithm> fb(max_thrs);
                        vector<ModelParameterVector> local_param_vec;
                        for (auto i = 0; i < max_thrs; ++i) { local_param_vec.push_back(new_param_vec); }

                        int tid = omp_get_thread_num();

                        size_t start_i = blocks[tid*2],
                               end_i = blocks[tid*2 + 1];

                        for (size_t i = start_i; i < end_i; ++i) {
                            if (good_clonotypes[i]) {
                                this->updateTempVec(fb[tid], maag_rep[i], local_param_vec[tid], error_mode);
                            }
                        }

                        for (size_t i = 0; i < new_param_vec.size(); ++i) {
                            #pragma omp atomic
                            new_param_vec[i] += local_param_vec[tid][i];
                        }
                    }
#else
                    MAAGForwardBackwardAlgorithm fb;
                    tp1 = std::chrono::system_clock::now();
                    for (size_t i = 0; i < rep_nonc.size(); ++i) {
                        if ((i+1) % 25000 == 0) {
                            cout << "Processed " << (int) i << " / " << (int) rep_nonc.size() << " MAAGs.\t" << endl;
                        }
                        if (good_clonotypes[i]) {
                            if(!this->updateTempVec(fb, maag_rep[i], new_param_vec, error_mode)) {
                                cout << "bad maag:\t" << (int) i << endl;
                            }
                        }
                    }
#endif
                }

                this->updateModel(rep_nonc, model, new_param_vec, maag_rep, prob_vec, prev_ll, removed, error_mode, memory_safe);

                std::cout << new_param_vec.error_prob() << std::endl;

                logLvec.push_back(prev_ll);
            }

            cout << endl << "Done. Resulting loglikelihood:\t" << loglikelihood(prob_vec) << endl << endl;

            return logLvec;
        }


    protected:

        bool updateTempVec(MAAGForwardBackwardAlgorithm &fb,
                           const MAAGnuc &maag,
                           ModelParameterVector &new_param_vec,
                           ErrorMode error_mode) const
        {
            if (!fb.process(maag, error_mode)) {
                return false;
            }

            while (!fb.is_empty()) {
                event_pair_t ep = fb.nextEvent();
                new_param_vec[ep.first] += ep.second;

                if (std::isnan(ep.second)) {
                    cout << "NaNs in the forw-back!" << endl;
                    return false;
                }
            }

            if (error_mode) { new_param_vec.set_error_prob(new_param_vec.error_prob() + fb.err_prob()); }

            if (maag.is_vj()) {
                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)] += fb.VJ_nuc_probs()[0];
                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 1)] += fb.VJ_nuc_probs()[1];
                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 2)] += fb.VJ_nuc_probs()[2];
                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 3)] += fb.VJ_nuc_probs()[3];
            } else {
                int k = 0;
                for (auto prev_nuc = 0; prev_nuc < 4; ++prev_nuc) {
                    for (auto next_nuc = 0; next_nuc < 4; ++next_nuc, ++k) {
                        new_param_vec[new_param_vec.event_index(VDJ_VAR_DIV_INS_NUC, prev_nuc, next_nuc)] += fb.VD_nuc_probs()[k];
                    }
                }

                k = 0;
                for (auto prev_nuc = 0; prev_nuc < 4; ++prev_nuc) {
                    for (auto next_nuc = 0; next_nuc < 4; ++next_nuc, ++k) {
                        new_param_vec[new_param_vec.event_index(VDJ_DIV_JOI_INS_NUC, prev_nuc, next_nuc)] += fb.DJ_nuc_probs()[k];
                    }
                }
            }

            return true;
        }


        void updateModel(const ClonesetViewNuc &rep_nonc,
                         ProbabilisticAssemblingModel &model,
                         ModelParameterVector &new_param_vec,
                         MAAGNucRepertoire &maag_rep,
                         vector<prob_t> &prob_vec,
                         prob_t &prev_ll,
                         size_t removed,
                         ErrorMode error_mode,
                         bool memory_safe) const
        {
            if (error_mode) { new_param_vec.set_error_prob(new_param_vec.error_prob() / (maag_rep.size() - removed)); }
//            new_param_vec.set_error_prob(.0003);

            new_param_vec.normaliseEventFamilies();

            model.updateModelParameterVector(new_param_vec);

            if (!memory_safe) {
                model.updateEventProbabilities(&maag_rep);
            }

            if (memory_safe) {
                prob_vec = model.computeFullProbabilities(rep_nonc, error_mode);
            } else {
                std::cout << "Recomputing generation probabilities...";
#ifdef USE_OMP
#pragma omp parallel for
#endif
                for (size_t i = 0; i < maag_rep.size(); ++i) {
                    prob_vec[i] = maag_rep[i].fullProbability();
                }
            }

            std::cout << " Done." << std::endl;

            prob_summary(prob_vec, prev_ll);
            prev_ll = loglikelihood(prob_vec);
        }

    };

}

#endif //YMIR_EM_ALGORITHM_H
