//
// Created by Vadim N. on 18/03/2015.
//

#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_


#include <chrono>
#include <ctime>

#include "Inference"


using namespace ymir;


#define YMIR_BENCHMARK(descr, expr) { tp1 = std::chrono::system_clock::now(); expr; tp2 = std::chrono::system_clock::now(); timepoints.emplace_back(descr, std::chrono::system_clock::to_time_t(tp2) - std::chrono::system_clock::to_time_t(tp1)); }


void write_vec(std::string filename, const std::vector<prob_t> &logLvec) {
    std::ofstream stream(filename);
    for (size_t i = 0; i < logLvec.size(); ++i) {
        stream << logLvec[i] << endl;
    }
}


//void infer_vdj_em(string descr, int niter, int sample, ErrorMode error_mode) {
//    YMIR_BENCHMARK(descr,
//                   logLvec = EMAlgorithm().statisticalInference(cloneset_vdj, vdj_single_model,
//                                                      EMAlgorithm::AlgorithmParameters()
//                                                              .set("niter", niter)
//                                                              .set("sample", sample),
//                                                      error_mode))
//    write_vec("/Users/vdn/Projects/ymir/benchmark/log/vdj_em_niter_" + to_string(niter)
//              + "_sample_" + to_string(sample)
//              + "_err_" + to_string(error_mode) +  ".txt", logLvec);
//}


int main(int argc, char* argv[]) {
    std::chrono::system_clock::time_point tp1, tp2;

    std::vector< std::pair<std::string, size_t> > timepoints;

    std::vector<prob_t> logLvec;

    time_t vj_single_parse, vdj_single_parse,
            vj_single_prob, vj_single_meta,
            vj_single_infer, vdj_single_prob,
            vdj_single_meta, vdj_single_infer;

    time_t vj_good_parse, vdj_good_parse,
            vj_good_prob, vj_good_meta,
            vj_good_infer, vdj_good_prob,
            vdj_good_meta, vdj_good_infer;
    
    std::string BENCH_DATA_FOLDER = argv[1];


    CDR3NucParser parser;

    string input_alpha_file = "alpha.250k.txt";
    string input_beta_file = "beta.250k.txt";


    //
    // TCR alpha chain repertoire - VJ recombination
    //
    VDJRecombinationGenes vj_single_genes("Vgene",
                                   BENCH_DATA_FOLDER + "trav.txt",
                                   "Jgene",
                                   BENCH_DATA_FOLDER + "traj.txt");

    Cloneset cloneset_vj;
    YMIR_BENCHMARK("Parsing VJ",
                   parser.openAndParse(BENCH_DATA_FOLDER + input_alpha_file,
                                       &cloneset_vj,
                                       vj_single_genes,
                                       NUCLEOTIDE,
                                       VJ_RECOMB,
                                       AlignmentColumnOptions(AlignmentColumnOptions::USE_PROVIDED, AlignmentColumnOptions::USE_PROVIDED),
                                       VDJAlignerParameters(2)))

    //
    // TCR beta chain repertoire - VDJ recombination
    //
    VDJRecombinationGenes vdj_single_genes("Vgene",
                                    BENCH_DATA_FOLDER + "trbv.txt",
                                    "Jgene",
                                    BENCH_DATA_FOLDER + "trbj.txt",
                                    "Dgene",
                                    BENCH_DATA_FOLDER + "trbd.txt");

    Cloneset cloneset_vdj;
    YMIR_BENCHMARK("Parsing VDJ",
                   parser.openAndParse(BENCH_DATA_FOLDER + input_beta_file,
                                       &cloneset_vdj,
                                       vdj_single_genes,
                                       NUCLEOTIDE,
                                       VDJ_RECOMB,
                                       AlignmentColumnOptions()
                                               .setV(AlignmentColumnOptions::USE_PROVIDED)
                                               .setD(AlignmentColumnOptions::OVERWRITE)
                                               .setJ(AlignmentColumnOptions::USE_PROVIDED),
                                       VDJAlignerParameters(3)))

    //
    // VJ MAAG
    //
    ProbabilisticAssemblingModel vj_single_model(BENCH_DATA_FOLDER + "../../models/hTRA"); //, EMPTY);

//    YMIR_BENCHMARK("VJ meta", vj_single_model.buildGraphs(cloneset_vj, SAVE_METADATA, NO_ERRORS))
//    YMIR_BENCHMARK("VJ prob", vj_single_model.computeFullProbabilities(cloneset_vj, NO_ERRORS, NUCLEOTIDE))


    //
    // VDJ MAAG
    //
    ProbabilisticAssemblingModel vdj_single_model(BENCH_DATA_FOLDER + "../../models/hTRB", EMPTY);

//    YMIR_BENCHMARK("VDJ meta", vdj_single_model.buildGraphs(cloneset_vdj, SAVE_METADATA, NO_ERRORS))
//    YMIR_BENCHMARK("VDJ prob", vdj_single_model.computeFullProbabilities(cloneset_vdj, NO_ERRORS, NUCLEOTIDE))


    //
    // VJ inference
    //
    // 10 - -4005118.35
//    EMAlgorithm().statisticalInference(cloneset_vj, vj_single_model);
//                                       EMAlgorithm::AlgorithmParameters().set("niter", 10),
//                                       NO_ERRORS);
    // 10 - -4036068.19
    //
//    SGAlgorithm().statisticalInference(cloneset_vj, vj_single_model);
//                                       EMAlgorithm::AlgorithmParameters().set("niter", 20));
//                                       NO_ERRORS);


    //
    // VDJ inference
    //
    int niter, sample, block;
    double alpha;
    niter = 4;
    sample = 50000;
    ErrorMode error_mode = NO_ERRORS;
//    YMIR_BENCHMARK("VDJ EM",
//                   logLvec = EMAlgorithm().statisticalInference(cloneset_vdj, vdj_single_model,
//                                                      EMAlgorithm::AlgorithmParameters()
//                                                              .set("niter", niter)
//                                                              .set("sample", sample),
//                                                      error_mode))
    write_vec("/Users/vdn/Projects/ymir/benchmark/log/vdj_em_niter_" + to_string(niter)
              + "_sample_" + to_string(sample)
              + "_err_" + to_string(error_mode) +  ".txt", logLvec);

    niter = 1;
    block = 5000;
    alpha = .6;
    sample = 50000;
    YMIR_BENCHMARK("VDJ SG",
                   logLvec = SGAlgorithm().statisticalInference(cloneset_vdj, vdj_single_model,
                                                      SGAlgorithm::AlgorithmParameters()
                                                              .set("niter", niter)
                                                              .set("block.size", block)
                                                              .set("alpha", alpha)
                                                              .set("prebuild", false)
                                                              .set("sample", sample),
                                                                error_mode))
    write_vec("/Users/vdn/Projects/ymir/benchmark/log/vdj_sg_niter_" + to_string(niter)
              + "_block_" + to_string(block)
              + "_alpha_" + to_string(alpha)
              + "_sample_" + to_string(sample)
              + "_err_" + to_string(error_mode) + ".txt", logLvec);


    //
    // Results
    //
    cout << "========================" << endl << "Results:" << endl;

    for (size_t i = 0; i < timepoints.size(); ++i) {
        cout << timepoints[i].first << ":\t" << timepoints[i].second << endl;
    }

    cout << endl << "========================" << endl;


    return 0;
}

#endif //_BENCHMARK_H_
