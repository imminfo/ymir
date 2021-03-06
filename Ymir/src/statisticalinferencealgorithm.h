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
                if (_json.value(param_name, "__NA__") == "__NA__") {
                    cout << "Obligatory parameter '" << param_name << "' hasn't been found, please re-run the algorithm with the supplied parameter." << endl;
                    return false;
                }
                return true;
            }


            json_value get(const string& param_name, const json_value& default_value) const {
                return _json.value(param_name, default_value);
            }


            AlgorithmParameters& set(const string& param_name, const json_value& value) {
                _json[param_name] = value;
                return *this;
            }


            const json_value& operator[](const string& param_name) const { return _json[param_name]; }


        private:

            json_value _json;

        };


        virtual ~StatisticalInferenceAlgorithm() { }

        virtual std::vector<prob_t> statisticalInference(const ClonesetViewNuc &repertoire,
                                                         ProbabilisticAssemblingModel &model,
                                                         const AlgorithmParameters &algo_param = AlgorithmParameters(),
                                                         ErrorMode error_mode = NO_ERRORS) const = 0;


        void filterOut(const ClonesetViewNuc &rep_nonc,
                       ProbabilisticAssemblingModel &model,
                       const MAAGNucRepertoire &maag_rep,
                       vector<prob_t> &prob_vec,
                       vector<bool> &good_clonotypes,
                       size_t &removed,
                       size_t &zero_prob,
                       size_t &no_alignments,
                       ErrorMode error_mode,
                       bool memory_safe) const
        {
            good_clonotypes.resize(rep_nonc.size(), true);
            prob_vec.resize(rep_nonc.size(), 0);
            removed = 0;
            zero_prob = 0;
            no_alignments = 0;

            if (!memory_safe) {
                cout << "Computing full assembling probabilities..." << endl;
#ifdef USE_OMP
#pragma omp parallel for
#endif
                for (size_t i = 0; i < maag_rep.size(); ++i) {
                    prob_vec[i] = maag_rep[i].fullProbability();
                }
            } else {
                prob_vec = model.computeFullProbabilities(rep_nonc, error_mode);
            }

            for (size_t i = 0; i < rep_nonc.size(); ++i) {
                if (rep_nonc[i].is_good()) {
                    if (std::isnan(prob_vec[i]) || (std::abs(prob_vec[i]) < 1e-80) || (std::abs(prob_vec[i]) >= 1)) {
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
