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
    class EMAlgorithm;
    class BlockEMAlgorithm;
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


            bool check(const string& param_name) const {
                return _json.get(param_name, "__NA__").asString() != "__NA__";
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
                cout << "\tBad alignments:\t" << (int) no_alignments << std::endl;
            } else {
                cout << "No clonotypes with error probabilities has been found. It's good in case you don't know." << std::endl;
            }
        }

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
                                          const AlgorithmParameters &algo_param = AlgorithmParameters().set("niter", 10),
                                          ErrorMode error_mode = NO_ERRORS) const {
            cout << "Statistical inference on a PAM:\t" << model.name() << endl;
            cout << "\tMurugan EM-algorithm.";
            if (error_mode == COMPUTE_ERRORS) {
                cout << "\t(with sequence errors)";
            }
            std::cout << std::endl;

            if (!algo_param.check("niter")) {
                cout << "Obligatory parameter 'niter' hasn't been found, please re-run the algorithm with the supplied parameter." << endl;
                return false;
            }

            ClonesetView rep_nonc = repertoire.noncoding().shuffle();
            cout << "Number of noncoding clonotypes:\t" << (size_t) rep_nonc.size() << endl;

            ModelParameterVector new_param_vec = model.event_probabilities();
            new_param_vec.fill(1);
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

                for (size_t i = 0; i < maag_rep.size(); ++i) {
                    if (i % 25000 == 0) {
                        cout << "Processed " << (int) i << " / " << (int) maag_rep.size() << " MAAGs.\t" << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << " sec." << endl;
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


    /**
    * \class BlockEMAlgorithm
    *
    * \brief Implementation of the online EM-algorithm for statistical inference of assembling model parameters.
    */
    class BlockEMAlgorithm : public EMAlgorithm {
    public:

        virtual bool statisticalInference(const ClonesetView& repertoire,
                                          ProbabilisticAssemblingModel & model,
                                          const AlgorithmParameters& algo_param = AlgorithmParameters().set("niter", 10).set("block.size", 2000),
                                          ErrorMode error_mode = NO_ERRORS) const
        {
            // shuffle input data at each step
            // subvec

            vector<size_t> indices;
            size_t start_i;
            size_t block_size = algo_param["block.size"].asUInt();
//            prob_t step = algo_param["block.size"].asDouble();
            ModelParameterVector new_param_vec = model.event_probabilities();
            vector<bool> changed(new_param_vec.size(), false);

            ClonesetView rep_nonc = repertoire.noncoding().shuffle();
            auto maag_rep = model.buildGraphs(rep_nonc, SAVE_METADATA, error_mode, NUCLEOTIDE, true);

            vector<prob_t> prob_vec(maag_rep.size(), 0);
            prob_t prev_ll = 0;

            MAAGForwardBackwardAlgorithm fb;
            for (size_t iter = 1; iter <= algo_param["niter"].asUInt(); ++iter) {

                new_param_vec.fill(0);

                if (start_i + block_size > maag_rep.size()) {
                    start_i = 0;
                    this->getRandomIndices(indices, maag_rep.size());
                }

                for (size_t maag_i = start_i; maag_i < std::min(maag_rep.size(), start_i + block_size); ++maag_i) {
                    // compute marginal probabilities for this block
                    // and update the temporary model parameter vector
                    this->updateTempVec(fb, maag_rep[indices[maag_i]], new_param_vec, changed, error_mode);
                }

                this->updateModel(model, new_param_vec, maag_rep, prob_vec, prev_ll, changed, error_mode);

                // update starting index
                start_i += block_size;
            }
        }

    protected:

        vector<size_t> getRandomIndices(vector<size_t> &indices, size_t size) const {

        }


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
                         ErrorMode error_mode) const
        {
            if (error_mode) { new_param_vec.set_error_prob(new_param_vec.error_prob() / (maag_rep.size())); }

            new_param_vec.normaliseEventFamilies();

            for (size_t i = 0; i < new_param_vec.size(); ++i) {
                if (!changed[i]) {
                    new_param_vec[i] = model.event_probabilities()[i];
                }
            }

            model.updateModelParameterVector(new_param_vec);
            model.updateEventProbabilities(&maag_rep);

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

#endif
