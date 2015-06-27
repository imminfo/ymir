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
    class StatisticalInferenceAlgorithm : protected ProbabilisticAssemblingModel {
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

        virtual bool statisticalInference(const ClonesetView &repertoire,
                                          ProbabilisticAssemblingModel &model,
                                          const AlgorithmParameters &algo_param = AlgorithmParameters().set("niter", 10)) const {

            cout << "Statistical inference on a PAM:\t" << model.name() << endl;

            ClonesetView rep_nonc = repertoire.noncoding();
            cout << "Number of noncoding clonotypes:\t" << (size_t) rep_nonc.size() << endl;

            MAAGRepertoire maag_rep = model.buildGraphs(rep_nonc);
            ModelParameterVector new_param_vec = model.event_probabilities();
            vector<prob_t> prob_vec;
            prob_vec.resize(rep_nonc.size(), 0);

            for (size_t iter = 1; iter <= algo_param["niter"].asUInt(); ++iter) {
                cout << "Iteration:\t" << (size_t) iter << endl;

                new_param_vec.clear();
                for (size_t i = 0; i < maag_rep.size(); ++i) {
                    MAAGForwardBackwardAlgorithm fb(maag_rep[i]);
                    while (!fb.is_empty()) {
                        event_pair_t ep = fb.nextEvent();
                        new_param_vec[ep.first] += ep.second;
                    }
                    prob_vec[i] = fb.fullProbability();
                }

                new_param_vec.normaliseEventFamilies();
                *(model._param_vec) = new_param_vec;
                model.make_builder();
                model.updateEventProbabilities(&maag_rep);

                cout << "Loglikelihood:\t" << loglikelihood(prob_vec) << endl;
            }

            model.make_assembler();
            return true;
        }

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