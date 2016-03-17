//
// Created by Vadim N. on 16/03/2016.
//

#ifndef YMIR_BENCHUTILS_H
#define YMIR_BENCHUTILS_H


#include <chrono>
#include <ctime>

#include "Inference"


using namespace ymir;
using namespace std;


#define YMIR_BENCHMARK(descr, expr) { tp1 = std::chrono::system_clock::now(); expr; tp2 = std::chrono::system_clock::now(); timepoints.emplace_back(descr, std::chrono::system_clock::to_time_t(tp2) - std::chrono::system_clock::to_time_t(tp1)); }


#define RUN_EM_INFERENCE(descr, cloneset, model, niter, sample, error_mode) { YMIR_BENCHMARK(descr, \
logLvec = EMAlgorithm().statisticalInference(cloneset, model, \
                                             EMAlgorithm::AlgorithmParameters() \
                                                     .set("niter", niter) \
                                                     .set("sample", sample), \
                                             error_mode)) \
temp_str  = std::string("/Users/vdn/Projects/ymir/benchmark/log/") + descr + "_em_niter_" + to_string(niter) \
+ "_sample_" + to_string(sample) \
+ "_err_" + to_string(error_mode); \
write_vec(temp_str + ".ll.txt", logLvec); \
write_vec(temp_str + ".pars.txt", model); \
  }


#define RUN_SG_INFERENCE(descr, cloneset, model, niter, block, alpha, beta, K, sample, error_mode) { YMIR_BENCHMARK(descr, \
logLvec = SGAlgorithm().statisticalInference(cloneset, model, \
                                             SGAlgorithm::AlgorithmParameters() \
                                                     .set("niter", niter) \
                                                     .set("block.size", block) \
                                                     .set("alpha", alpha) \
                                                     .set("beta", beta) \
                                                     .set("K", K) \
                                                     .set("prebuild", false) \
                                                     .set("sample", sample), \
                                             error_mode)) \
temp_str = std::string("/Users/vdn/Projects/ymir/benchmark/log/") + descr + "_sg_niter_" + to_string(niter) \
+ "_block_" + to_string(block) \
+ "_alpha_" + to_string(alpha) \
+ "_beta_" + to_string(beta) \
+ "_K_" + to_string(K) \
+ "_sample_" + to_string(sample) \
+ "_err_" + to_string(error_mode); \
write_vec(temp_str + ".ll.txt", logLvec); \
write_vec(temp_str + ".pars.txt", model); \
  }


void write_vec(std::string filename, const std::vector<prob_t> &logLvec) {
    std::ofstream stream(filename);
    for (size_t i = 0; i < logLvec.size(); ++i) {
        stream << logLvec[i] << endl;
    }
}


void write_vec(std::string filename, const ProbabilisticAssemblingModel &model) {
    std::ofstream stream(filename);
    for (size_t i = 0; i < model.event_probabilities().size(); ++i) {
        stream << model.event_probabilities()[i] << endl;
    }
    stream << model.event_probabilities().error_prob() << endl;
}


#endif //YMIR_BENCHUTILS_H
