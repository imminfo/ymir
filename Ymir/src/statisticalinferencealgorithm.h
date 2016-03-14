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


#include <chrono>
#include <ctime>

#include "probabilisticassemblingmodel.h"
#include "maagforwardbackwardalgorithm.h"
#include "tools.h"


namespace ymir {

    class StatisticalInferenceAlgorithm;


    /**
    * \class StatisticalInferenceAlgorithm
    *
    * \brief Interface for algorithms for statistical inference of assembling model parameters.
    */
    class StatisticalInferenceAlgorithm {
    public:

        struct AlgorithmParameters {


            AlgorithmParameters() {}


            bool check(const string& param_name) const {
                if (_json.get(param_name, "__NA__").asString() == "__NA__") {
                    cout << "Obligatory parameter '" << param_name << "' hasn't been found, please re-run the algorithm with the supplied parameter." << endl;
                    return false;
                }
                return true;
            }


            AlgorithmParameters& set(const string& param_name, const Json::Value& value) {
                _json[param_name] = value;
                return *this;
            }


            const Json::Value& operator[](const string& param_name) const { return _json[param_name]; }


        private:

            Json::Value _json;

        };


        virtual ~StatisticalInferenceAlgorithm() { }

        virtual bool statisticalInference(const ClonesetView &repertoire,
                                          ProbabilisticAssemblingModel &model,
                                          const AlgorithmParameters &algo_param = AlgorithmParameters(),
                                          ErrorMode error_mode = NO_ERRORS) const =0;


        void filterOut(const ClonesetView &rep_nonc,
                       const MAAGRepertoire &maag_rep,
                       vector<prob_t> &prob_vec,
                       vector<bool> &good_clonotypes,
                       size_t &removed,
                       size_t &zero_prob,
                       size_t &no_alignments) const
        {
            good_clonotypes.resize(maag_rep.size(), true);
            prob_vec.resize(maag_rep.size(), 0);
            removed = 0;
            zero_prob = 0;
            no_alignments = 0;

            for (size_t i = 0; i < maag_rep.size(); ++i) {
                if (rep_nonc[i].is_good()) {
                    prob_vec[i] = maag_rep[i].fullProbability();
                    if (std::isnan(prob_vec[i]) || prob_vec[i] == 0) {
                        good_clonotypes[i] = false;
                        ++removed;
                        ++zero_prob;
                    }
                } else {
                    good_clonotypes[i] = false;
                    ++removed;
                    ++no_alignments;
                }
            }

            if (removed) {
                cout << "Removed " << (int) removed
                << " error-probability clonotypes. Check your minimal Diversity gene length to align and other parameters to make sure it won't happen again in the future." << endl;
                cout << "\tZero probabilities:\t" << (int) zero_prob << std::endl;
                cout << "\tBad alignments:\t" << (int) no_alignments << " (replaces with zeros)" << std::endl;
            } else {
                cout << "No clonotypes with error probabilities has been found. It's good in case you don't know." << std::endl;
            }
        }

    };

}

#endif
