/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vadim dot nazarov at mailbox dot com>
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


#include "assemblingstatisticalmodel.h"

#include "unordered_map"


namespace ymir {

    class StatisticalInferenceAlgorithm;
    class EMAlgorithm;
    class OnlineEMAlgorithm;


    /**
    * \class StatisticalInferenceAlgorithm
    *
    * \brief Interface for algorithms for statistical inference of assembling model parameters.
    */
    class StatisticalInferenceAlgorithm : protected AssemblingStatisticalModel {
    public:

        struct AlgorithmParameters {


            AlgorithmParameters() {}


            AlgorithmParameters& set(const string& param_name, const Json::Value& value) {
                _json[param_name] = value;
                return *this;
            }


            Json::Value& operator[](const string& param_name) { return _json[param_name]; }


        private:

            Json::Value _json;

        };


        virtual bool statisticalInference(const ClonalRepertoireView& repertoire,
                AssemblingStatisticalModel& model,
                const AlgorithmParameters& algo_param = AlgorithmParameters()) const =0;


        // function for aggregating P list(s) / ModelParameterVector from assembly scenario matrices

        // function for check if all sequences in clonal repertoire is nucleotide
    };


    /**
    * \class EMAlgorithm
    *
    * \brief Implementation of EM-algorithm for statistical inference of assembling model parameters.
    */
    class EMAlgorithm : public StatisticalInferenceAlgorithm {
    public:

        virtual bool statisticalInference(const ClonalRepertoireView& repertoire,
                AssemblingStatisticalModel& model,
                const AlgorithmParameters& algo_param = AlgorithmParameters().set("niter", 10)) const {

        }
    };


    /**
    * \class OnlineEMAlgorithm
    *
    * \brief Implementation of online EM-algorithm for statistical inference of assembling model parameters.
    */
    class OnlineEMAlgorithm : public StatisticalInferenceAlgorithm {
    public:

        virtual bool statisticalInference(const ClonalRepertoireView& repertoire,
                AssemblingStatisticalModel& model,
                const AlgorithmParameters& algo_param = AlgorithmParameters().set("niter", 10).set("step.size", 2000)) const {

        }
    };
}

#endif