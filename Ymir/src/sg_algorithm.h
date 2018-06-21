// //
// // Created by Vadim N. on 11/03/2016.
// //

// #ifndef YMIR_SG_ALGORITHM_H
// #define YMIR_SG_ALGORITHM_H


#include "statisticalinferencealgorithm.h"


// namespace ymir {

//     class SGAlgorithm;


//     /**
//    * \class BlockEMAlgorithm
//    *
//    * \brief Implementation of the online EM-algorithm for statistical inference of assembling model parameters.
//    */
//    class SGAlgorithm : public StatisticalInferenceAlgorithm {
//    public:

//        virtual std::vector<prob_t> statisticalInference(const ClonesetViewNuc& repertoire,
//                                                         ProbabilisticAssemblingModel & model,
//                                                         const AlgorithmParameters& algo_param = AlgorithmParameters()
//                                                                 .set("niter", 10)
//                                                                 .set("memory-safe", false)
//                                                                 .set("niter", 50)
//                                                                 .set("block.size", 5000)
//                                                                 .set("parallel.blocks", 2)
//                                                                 .set("alpha", .6)
//                                                                 .set("beta", 1.)
//                                                                 .set("K", 2.),
//                                                         ErrorMode error_mode = NO_ERRORS) const
//        {
//            cout << "Statistical inference on a PAM:\t" << model.name() << endl;
//            cout << "\tOnline EM-algorithm.";
//            if (error_mode == COMPUTE_ERRORS) {
//                cout << "\t(with sequence errors)";
//            }
//            std::cout << std::endl;


//            if (!algo_param.check("niter")
//                && !algo_param.check("memory-safe")
//                && !algo_param.check("block.size")
//                && !algo_param.check("parallel.blocks")
//                && !algo_param.check("alpha")
//                && !algo_param.check("beta")
//                && !algo_param.check("K"))
//            {
//                return std::vector<prob_t>();
//            }

//            std::cout << "\t -- #iterations: " << (size_t) algo_param["niter"].asUInt() << std::endl;
//            std::cout << "\t -- block size:  " << (size_t) algo_param["block.size"].asUInt() << std::endl;
//            std::cout << "\t -- par. blocks: " << (size_t) algo_param["parallel.blocks"].asUInt() << std::endl;
//            std::cout << "\t -- alpha:       " << (double) algo_param["alpha"].asDouble() << std::endl;
//            std::cout << "\t -- beta:        " << (double) algo_param["beta"].asDouble() << std::endl;
//            std::cout << "\t -- K:           " << (double) algo_param["K"].asDouble() << std::endl;

//            bool memory_safe = algo_param["memory-safe"].asBool();
//            if (memory_safe) {
//                std::cout << "\tworking in the memory-safe mode: rebuild graphs at each step"
//                          << std::endl;
//            } else {
//                std::cout << "\tworking in the speed-wise mode: build graphs only once"
//                          << std::endl;
//            }

//            // size_t sample = algo_param["sample"].asUInt();
//            ClonesetViewNuc rep_nonc = repertoire.noncoding();
//            cout << "Number of noncoding clonotypes:\t" << (size_t) rep_nonc.size() << endl;

//            std::vector<prob_t> logLvec;


//            size_t start_i = rep_nonc.size();
//            size_t block_size = algo_param["block.size"].asUInt();
//            size_t blocks_in_parallel = algo_param["parallel.blocks"].asUInt();
//            prob_t alpha = algo_param["alpha"].asDouble(); // step(k) = (k + 2)^(-alpha), .5 < alpha <= 1
//            prob_t beta = algo_param["beta"].asDouble();
//            prob_t Kparam = algo_param["K"].asDouble(); // step(k) = (k + 2)^(-alpha), .5 < alpha <= 1

//            ModelParameterVector new_param_vec = model.event_probabilities();
//            new_param_vec.fill(1);
//            // Error probability
//            new_param_vec.set_error_prob(.0003);
//            new_param_vec.normaliseEventFamilies();
//            model.updateModelParameterVector(new_param_vec);

//            std::vector<bool> changed(new_param_vec.size(), false);

//            MAAGNucRepertoire maag_rep;
//            if (!memory_safe) {
//                 auto maag_rep = model.buildGraphs(rep_nonc, SAVE_METADATA, error_mode, true);
//            }

//            vector<prob_t> prob_vec;

//             prob_t prev_ll = 0, cur_ll = 0;

//             vector<bool> good_clonotypes;
//             size_t removed, zero_prob, no_alignments;
//             this->filterOut(rep_nonc, model, maag_rep, prob_vec, good_clonotypes, removed, zero_prob, no_alignments, error_mode, memory_safe);

//             cout << endl << "Initial data summary:" << endl;
//             prob_summary(prob_vec);
//             std::cout << model.event_probabilities().error_prob() << std::endl;
//             prev_ll = loglikelihood(prob_vec);
//             logLvec.push_back(prev_ll);



//             std::cout << "Infer parameters..." << std::endl;
// #ifdef USE_OMP
//             auto max_thrs = omp_get_max_threads();
//             // how many threads work per block
//             auto thread_per_block = max_thrs / blocks_in_parallel

//             std::vector<size_t> blocks;
//             size_t block_step = std::min(rep_nonc.size() / max_thrs + 1, rep_nonc.size());
//             blocks.push_back(0);
//             blocks.push_back(block_step);
//             for (auto i = 1; i < max_thrs; ++i) {
//                 blocks.push_back(*--blocks.end());
//                 blocks.push_back(std::min(*--blocks.end() + block_step, rep_nonc.size()));
//             }

//             #pragma omp parallel
//             {
//                 vector<MAAGForwardBackwardAlgorithm> fb(max_thrs);
//                 vector<ModelParameterVector> local_param_vec;
//                 for (auto i = 0; i < max_thrs; ++i) { local_param_vec.push_back(new_param_vec); }

//                 int tid = omp_get_thread_num();

//                 size_t start_i = blocks[tid*2],
//                        end_i = blocks[tid*2 + 1];

//                 for (size_t i = start_i; i < end_i; ++i) {
//                     if (good_clonotypes[i]) {
//                         this->updateTempVec(fb[tid], model.buildGraphs(rep_nonc[i], SAVE_METADATA, error_mode), local_param_vec[tid], error_mode);
//                     }
//                 }

//                 for (size_t i = 0; i < new_param_vec.size(); ++i) {
//                     #pragma omp atomic
//                     new_param_vec[i] += local_param_vec[tid][i];
//                 }
//             }
// #else
//             MAAGForwardBackwardAlgorithm fb;
//             tp1 = std::chrono::system_clock::now();
//             for (size_t i = 0; i < rep_nonc.size(); ++i) {
//                 if ((i+1) % 25000 == 0) {
//                     cout << "Processed " << (int) i << " / " << (int) rep_nonc.size() << " MAAGs.\t" << endl;
//                 }
//                 if (good_clonotypes[i]) {
//                     if(!this->updateTempVec(fb, model.buildGraphs(rep_nonc[i], SAVE_METADATA, error_mode), new_param_vec, error_mode)) {
//                         cout << "bad maag:\t" << (int) i << endl;
//                     }
//                 }
//             }
// #endif



//            MAAGForwardBackwardAlgorithm fb;
//            for (size_t iter = 1; iter <= algo_param["niter"].asUInt(); ++iter) {
//                if (start_i + block_size > indices.size()) {
//                    start_i = 0;
//                    std::random_shuffle(indices.begin(), indices.end());
//                } else {
//                    start_i += block_size;
//                }

//                std::cout << "=======================" << std::endl
//                << "Iteration: " << (size_t) iter << " Block: [" << (int) start_i << ":" << (int) std::min(indices.size() - 1, start_i + block_size - 1)
//                << "]" << std::endl << "=======================" << std::endl;

//                new_param_vec.fill(0);

//                for (size_t maag_i = start_i; maag_i < std::min(indices.size(), start_i + block_size); ++maag_i) {
//                    // compute marginal probabilities for this block
//                    // and update the temporary model parameter vector
//                    this->updateTempVec(fb, maag_rep[indices[maag_i]], new_param_vec, changed, error_mode);
//                }

//                std::cout << "Err:" << new_param_vec.error_prob() << std::endl;
//                this->updateModel(model, new_param_vec, maag_rep, prob_vec, prev_ll, changed, pow(beta*iter + Kparam, -alpha), error_mode);

//                logLvec.push_back(prev_ll);
//            }

//            return logLvec;
//        }

//    protected:

//        void updateTempVec(MAAGForwardBackwardAlgorithm &fb,
//                           const MAAGnuc &maag,
//                           ModelParameterVector &new_param_vec,
//                           vector<bool> &changed,
//                           ErrorMode error_mode) const
//        {
//            if (!fb.process(maag, error_mode)) {
//                 return false;
//             }

//             while (!fb.is_empty()) {
//                 event_pair_t ep = fb.nextEvent();
//                 new_param_vec[ep.first] += ep.second;

//                 if (std::isnan(ep.second)) {
//                     cout << "NaNs in the forw-back!" << endl;
//                     return false;
//                 }
//             }

//            if (error_mode) { new_param_vec.set_error_prob(new_param_vec.error_prob() + fb.err_prob()); }

//            if (maag.is_vj()) {
//                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)] += fb.VJ_nuc_probs()[0];
//                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 1)] += fb.VJ_nuc_probs()[1];
//                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 2)] += fb.VJ_nuc_probs()[2];
//                new_param_vec[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 3)] += fb.VJ_nuc_probs()[3];

//                changed[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)] = true;
//                changed[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 1)] = true;
//                changed[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 2)] = true;
//                changed[new_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 3)] = true;
//            } else {
//                int k = 0;
//                 for (auto prev_nuc = 0; prev_nuc < 4; ++prev_nuc) {
//                     for (auto next_nuc = 0; next_nuc < 4; ++next_nuc, ++k) {
//                         new_param_vec[new_param_vec.event_index(VDJ_VAR_DIV_INS_NUC, prev_nuc, next_nuc)] += fb.VD_nuc_probs()[k];
//                         changed[new_param_vec.event_index(VDJ_VAR_DIV_INS_NUC, prev_nuc, next_nuc)] = true;
//                     }
//                 }

//                 k = 0;
//                 for (auto prev_nuc = 0; prev_nuc < 4; ++prev_nuc) {
//                     for (auto next_nuc = 0; next_nuc < 4; ++next_nuc, ++k) {
//                         new_param_vec[new_param_vec.event_index(VDJ_DIV_JOI_INS_NUC, prev_nuc, next_nuc)] += fb.DJ_nuc_probs()[k];
//                         changed[new_param_vec.event_index(VDJ_DIV_JOI_INS_NUC, prev_nuc, next_nuc)] = true
//                     }
//                 }
//            }
//        }


//        void updateModel(const ClonesetViewNuc &rep_nonc,
//                         ProbabilisticAssemblingModel &model,
//                         vector<ModelParameterVector> &new_param_vec_storage,
//                         MAAGNucRepertoire &maag_rep,
//                         vector<prob_t> &prob_vec,
//                         vector<vector<bool> > &changed_storage,
//                         prob_t step_k,
//                         prob_t &prev_ll,
//                         size_t removed,
//                         ErrorMode error_mode) const
//        {
// //            if (error_mode) { new_param_vec.set_error_prob(new_param_vec.error_prob() / maag_rep.size()); }
//            if (error_mode) { new_param_vec.set_error_prob(.0003); }

//            vector<size_t> n_values_per_param(new_param_vec_storage[0].size());
//            auto final_param_vec = new_param_vec_storage[0];
//            final_param_vec.fill(0);

//            for (size_t storage_i = 0; storage_i < new_param_vec_storage.size(); ++storage_i) {
//               new_param_vec_storage[storage_i].normaliseEventFamilies();

//               for (size_t i = 0; i < new_param_vec_storage[storage_i].size(); ++i) {
//                 n_values_per_param[i] += static_cast<size_t>(changed_storage[storage_i][i]);

//                 if (!changed_storage[storage_i][i]) {
//                    final_param_vec[i] = model.event_probabilities()[i];
//                 } else {
//                    final_param_vec[i] = step_k * model.event_probabilities()[i] + (1 - step_k) * new_param_vec_storage[storage_i][i];
//                 }
//               }
//            }

//            for (size_t i = 0; i < final_param_vec.size(); ++i) {
//               // Take means from each changed values
//               if (n_values_per_param[i]) {
//                 final_param_vec[i] /= n_values_per_param[i];
//               }
//            }

//            final_param_vec.normaliseEventFamilies();

//            model.updateModelParameterVector(final_param_vec);

//            if (!memory_safe) {
//                 model.updateEventProbabilities(&maag_rep);
//             }

//            for (size_t i = 0; i < maag_rep.size(); ++i) {
//                prob_vec[i] = maag_rep[i].fullProbability();
//            }
//            prob_summary(prob_vec, prev_ll);
//            prev_ll = loglikelihood(prob_vec);

//            for (size_t storage_i = 0; storage_i < changed_storage.size(); ++storage_i) {
//                 changed_storage[storage_i].clear();
//                 changed_storage[storage_i].resize(final_param_vec.size(), false);
//            }
           
//        }

//    };

// }

// #endif //YMIR_SG_ALGORITHM_H
