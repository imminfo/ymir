//
// Created by Vadim N. on 11/03/2016.
//

#ifndef YMIR_EM_ALGORITHM_H
#define YMIR_EM_ALGORITHM_H


#include "statisticalinferencealgorithm.h"
#ifdef USE_OMP
#include <omp.h>
#endif


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
                                                         const AlgorithmParameters &algo_param = AlgorithmParameters().set("niter", 10).set("sample", 50000),
                                                         ErrorMode error_mode = NO_ERRORS) const {
            cout << "Statistical inference on a PAM:\t" << model.name() << endl;
            cout << "\tMurugan EM-algorithm.";
            if (error_mode == COMPUTE_ERRORS) {
                cout << "\t(with sequence errors)";
            }
            std::cout << std::endl;

            if (!algo_param.check("niter") && !algo_param.check("sample")) {
                return std::vector<prob_t>();
            }

            cout << "\t# iterations"
                 << (int) algo_param["niter"].asUInt()
                 << std::endl;

            size_t sample = algo_param["sample"].asUInt();
            if (sample) {
                std::cout << "\tsample size (doesn't work yet): "
                          << (int) sample
                          << std::endl;
            }
            ClonesetViewNuc rep_nonc = repertoire.noncoding();
            cout << "Number of noncoding clonotypes:\t" << (size_t) rep_nonc.size() << endl;

            std::vector<prob_t> logLvec;


            ModelParameterVector new_param_vec = model.event_probabilities();
            new_param_vec.fill(1);


            //
            // SET ERROR HERE
            //
            new_param_vec.set_error_prob(.0003);


            new_param_vec.normaliseEventFamilies();
            model.updateModelParameterVector(new_param_vec);

            MAAGNucRepertoire maag_rep = model.buildGraphs(rep_nonc, SAVE_METADATA, error_mode, true);

            vector<prob_t> prob_vec;

            prob_t prev_ll = 0, cur_ll = 0;

            cout << "Computing full assembling probabilities..." << endl;
            vector<bool> good_clonotypes;
            size_t removed, zero_prob, no_alignments;
            this->filterOut(rep_nonc, maag_rep, prob_vec, good_clonotypes, removed, zero_prob, no_alignments);

            cout << endl << "Initial data summary:" << endl;
            prob_summary(prob_vec);
            std::cout << model.event_probabilities().error_prob() << std::endl;
            prev_ll = loglikelihood(prob_vec);
            logLvec.push_back(prev_ll);

            std::chrono::system_clock::time_point tp1, tp2;
            MAAGForwardBackwardAlgorithm fb;
            for (size_t iter = 1; iter <= algo_param["niter"].asUInt(); ++iter) {
                cout << endl << "Iteration:\t" << (size_t) iter << " / " << (size_t) algo_param["niter"].asUInt() << endl;

                new_param_vec.fill(0);

#ifdef USE_OMP
#if OMP_THREADS == -1
            #pragma omp parallel
#else
            #pragma omp parallel num_threads(OMP_THREADS)
#endif
                {
                    auto local_param_vec = new_param_vec;
                    int tid = omp_get_thread_num();

                    size_t start_i = 0, end_i = 0;
                    #pragma omp for nowait
                    for (size_t i = start_i; i < end_i; ++i) {
                        this->updateTempVec(fb, maag_rep[i], local_param_vec, error_mode);
                    }

                    for (size_t i = 0; i < new_param_vec.size(); ++i) {
                        #pragma omp atomic
                        new_param_vec[i] += local_param_vec[i];
                    }
                }
#else

                tp1 = std::chrono::system_clock::now();
                for (size_t i = 0; i < maag_rep.size(); ++i) {
                    if ((i+1) % 25000 == 0) {
                        cout << "Processed " << (int) i << " / " << (int) maag_rep.size() << " MAAGs.\t" << (std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) - std::chrono::system_clock::to_time_t(tp1)) << " sec." << endl;
                    }
                    if (good_clonotypes[i]) {
                        if(!this->updateTempVec(fb, maag_rep[i], new_param_vec, error_mode)) {
                            cout << "bad maag:\t" << (int) i << endl;
                        }
                    }
                }
#endif

                this->updateModel(model, new_param_vec, maag_rep, prob_vec, prev_ll, removed, error_mode);

                std::cout << new_param_vec.error_prob() << std::endl;

                logLvec.push_back(prev_ll);
            }

            cout << endl << "Done. Resulting loglikelihood:\t" << loglikelihood(prob_vec) << endl << endl;

            return logLvec;
        }


    protected:

        bool updateTempVec(MAAGForwardBackwardAlgorithm &fb,
                           MAAGnuc &maag,
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


        void updateModel(ProbabilisticAssemblingModel &model,
                         ModelParameterVector &new_param_vec,
                         MAAGNucRepertoire &maag_rep,
                         vector<prob_t> &prob_vec,
                         prob_t &prev_ll,
                         size_t removed,
                         ErrorMode error_mode) const
        {
            if (error_mode) { new_param_vec.set_error_prob(new_param_vec.error_prob() / (maag_rep.size() - removed)); }
//            new_param_vec.set_error_prob(.0003);

            new_param_vec.normaliseEventFamilies();

            model.updateModelParameterVector(new_param_vec);
            model.updateEventProbabilities(&maag_rep);

            std::cout << "Recomputing generation probabilities...";
#ifdef USE_OMP
            #if OMP_THREADS == -1
            #pragma omp parallel for
#else
            #pragma omp parallel for num_threads(OMP_THREADS)
#endif
#endif
            for (size_t i = 0; i < maag_rep.size(); ++i) {
                prob_vec[i] = maag_rep[i].fullProbability();
            }
            std::cout << " Done." << std::endl;

            prob_summary(prob_vec, prev_ll);
            prev_ll = loglikelihood(prob_vec);
        }

    };

}

#endif //YMIR_EM_ALGORITHM_H
