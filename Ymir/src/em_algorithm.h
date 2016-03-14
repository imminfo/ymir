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
        virtual bool statisticalInference(const ClonesetView &repertoire,
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
                return false;
            }


            size_t sample = algo_param["sample"].asUInt();
            ClonesetView rep_nonc = repertoire.noncoding().sample(sample);
            cout << "Number of noncoding clonotypes:\t" << (size_t) rep_nonc.size() << endl;


            ModelParameterVector new_param_vec = model.event_probabilities();
            new_param_vec.fill(1);
//            new_param_vec.set_error_prob(.003);
            new_param_vec.normaliseEventFamilies();
            model.updateModelParameterVector(new_param_vec);

            MAAGRepertoire maag_rep = model.buildGraphs(rep_nonc, SAVE_METADATA, error_mode, NUCLEOTIDE, true);

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

            std::chrono::system_clock::time_point tp1, tp2;
            MAAGForwardBackwardAlgorithm fb;
            for (size_t iter = 1; iter <= algo_param["niter"].asUInt(); ++iter) {
                cout << endl << "Iteration:\t" << (size_t) iter << endl;

                new_param_vec.fill(0);

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

                this->updateModel(model, new_param_vec, maag_rep, prob_vec, prev_ll, removed, error_mode);

                std::cout << new_param_vec.error_prob() << std::endl;

            }

            cout << endl << "Done. Resulting loglikelihood:\t" << loglikelihood(prob_vec) << endl << endl;

            return true;
        }


    protected:

        bool updateTempVec(MAAGForwardBackwardAlgorithm &fb,
                           MAAG &maag,
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
                         MAAGRepertoire &maag_rep,
                         vector<prob_t> &prob_vec,
                         prob_t &prev_ll,
                         size_t removed,
                         ErrorMode error_mode) const
        {
            if (error_mode) { new_param_vec.set_error_prob(new_param_vec.error_prob() / (maag_rep.size() - removed)); }
            new_param_vec.normaliseEventFamilies();

            model.updateModelParameterVector(new_param_vec);
            model.updateEventProbabilities(&maag_rep);

            for (size_t i = 0; i < maag_rep.size(); ++i) {
                prob_vec[i] = maag_rep[i].fullProbability();
            }
            prob_summary(prob_vec, prev_ll);
            prev_ll = loglikelihood(prob_vec);
        }

    };

}

#endif //YMIR_EM_ALGORITHM_H
