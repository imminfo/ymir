/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vdn at mailbox dot com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _STATISTICALINFERENCEALGORITHM_H
#define _STATISTICALINFERENCEALGORITHM_H


#include "probabilisticassemblingmodel.h"
#include "maagforwardbackwardalgorithm.h"
#include "tools.h"


namespace ymir {

    class StatisticalInferenceAlgorithm;
    class EMAlgorithm;
    class OnlineEMAlgorithm;


    /**
    * \class StatisticalInferenceAlgorithm
    *
    * \brief Interface for algorithms for statistical inference of assembling model parameters.
    */
    class StatisticalInferenceAlgorithm {
    public:

        struct AlgorithmParameters {


            AlgorithmParameters() {}


            AlgorithmParameters& set(const string& param_name, const Json::Value& value) {
                _json[param_name] = value;
                return *this;
            }


            const Json::Value& operator[](const string& param_name) const { return _json[param_name]; }


        private:

            Json::Value _json;

        };


        virtual bool statisticalInference(const ClonesetView &repertoire,
                                          ProbabilisticAssemblingModel &model,
                                          const AlgorithmParameters &algo_param = AlgorithmParameters()) const =0;

    };


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
                                          const AlgorithmParameters &algo_param = AlgorithmParameters().set("niter", 10)) const {

            cout << "Statistical inference on a PAM:\t" << model.name() << endl;

            ClonesetView rep_nonc = repertoire.noncoding();
            cout << "Number of noncoding clonotypes:\t" << (size_t) rep_nonc.size() << endl;

            ModelParameterVector new_param_vec = model.event_probabilities();
            new_param_vec.fill(1);
            new_param_vec.normaliseEventFamilies();
            model.updateEventProbabilitiesVector(new_param_vec);

            cout << "Building MAAGs..." << endl;
            MAAGRepertoire maag_rep = model.buildGraphs(rep_nonc, SAVE_METADATA, false, false);
//            if (maag_rep[0].is_vj()) { maag_rep.erase(maag_rep.begin() + 61877); }

            vector<prob_t> prob_vec;
            prob_vec.resize(maag_rep.size(), 0);

            prob_t prev_ll = 0, cur_ll = 0;

            cout << endl << "Initial data summary:" << endl;
            for (size_t i = 0; i < maag_rep.size(); ++i) {
                prob_vec[i] = maag_rep[i].fullProbability();
            }
            prob_summary(prob_vec);
            prev_ll = loglikelihood(prob_vec);

            for (size_t iter = 1; iter <= algo_param["niter"].asUInt(); ++iter) {
                cout << endl << "Iteration:\t" << (size_t) iter << endl;

                new_param_vec.fill(0);

                for (size_t i = 0; i < maag_rep.size(); ++i) {
                    MAAGForwardBackwardAlgorithm fb(maag_rep[i]);

                    while (!fb.is_empty()) {
                        event_pair_t ep = fb.nextEvent();
                        new_param_vec[ep.first] += ep.second;
                    }

                    if (maag_rep[i].is_vj()) {
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
                }

                new_param_vec.normaliseEventFamilies();

                model.updateEventProbabilitiesVector(new_param_vec);
                model.updateEventProbabilities(&maag_rep);

                for (size_t i = 0; i < maag_rep.size(); ++i) {
                    prob_vec[i] = maag_rep[i].fullProbability();
                }
                prob_summary(prob_vec, prev_ll);
                prev_ll = loglikelihood(prob_vec);
            }

            return true;
        }


        /**
         *
         */
        void iteration() const
        {}

    };


    /**
    * \class OnlineEMAlgorithm
    *
    * \brief Implementation of the online EM-algorithm for statistical inference of assembling model parameters.
    */
    class OnlineEMAlgorithm : public StatisticalInferenceAlgorithm {
    public:

        virtual bool statisticalInference(const ClonesetView& repertoire,
                ProbabilisticAssemblingModel & model,
                const AlgorithmParameters& algo_param = AlgorithmParameters().set("niter", 10).set("step.size", 2000)) const {

        }
    };
}

#endif