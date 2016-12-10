//
// Created by Vadim N. on 18/03/2015.
//

#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_


#include "benchutils.h"


int main(int argc, char* argv[]) {

    if (argc < 2) {
      std::cout << "Please re-run the script with the data folder path supplied." << std::endl;
      return 1;
    }

    std::chrono::system_clock::time_point tp1, tp2;

    std::vector< std::pair<std::string, size_t> > timepoints;

    std::vector<prob_t> logLvec;

    std::string temp_str;
    
    std::string BENCH_DATA_FOLDER = argv[1];

    ClonesetNuc cloneset_vj, cloneset_vdj, cloneset_vj_noncoding, cloneset_vdj_noncoding;


    string input_alpha_file = "alpha.full.500k.txt", input_alpha_file_nonc = "alpha.noncoding.100k.txt";
    string input_beta_file = "beta.noncoding.500k.txt", input_beta_file_nonc = "beta.noncoding.100k.txt";


    auto vdj_aligner_parameters_nuc = VDJAlignerParameters(3,
                                                           VDJAlignmentEventScore(AlignmentEventScore(1, -1, 1),
                                                                                  AlignmentEventScore(1, -1, 1),
                                                                                  AlignmentEventScore(1, -1, 1)),
                                                           VDJAlignmentScoreThreshold(6, 3, 5));

    ParserNuc parser(new NaiveCDR3NucleotideAligner());


    //
    // TCR alpha chain repertoire - VJ recombination
    //
    VDJRecombinationGenes vj_single_genes("Vgene",
                                   BENCH_DATA_FOLDER + "trav.txt",
                                   "Jgene",
                                   BENCH_DATA_FOLDER + "traj.txt");

    YMIR_BENCHMARK("Parsing VJ",
                   parser.openAndParse(BENCH_DATA_FOLDER + input_alpha_file,
                                       &cloneset_vj,
                                       vj_single_genes,
                                       VJ_RECOMB,
                                       AlignmentColumnOptions(AlignmentColumnOptions::OVERWRITE,
                                                              AlignmentColumnOptions::OVERWRITE),
                                       vdj_aligner_parameters_nuc))

    //
    // TCR beta chain repertoire - VDJ recombination
    //
    VDJRecombinationGenes vdj_single_genes("Vgene",
                                    BENCH_DATA_FOLDER + "trbv.txt",
                                    "Jgene",
                                    BENCH_DATA_FOLDER + "trbj.txt",
                                    "Dgene",
                                    BENCH_DATA_FOLDER + "trbd.txt");

    YMIR_BENCHMARK("Parsing VDJ",
                   parser.openAndParse(BENCH_DATA_FOLDER + input_beta_file,
                                       &cloneset_vdj,
                                       vdj_single_genes,
                                       VDJ_RECOMB,
                                       AlignmentColumnOptions(AlignmentColumnOptions::OVERWRITE,
                                                              AlignmentColumnOptions::OVERWRITE,
                                                              AlignmentColumnOptions::OVERWRITE),
                                       vdj_aligner_parameters_nuc))

    //
    // VJ MAAG
    //
    vector<prob_t> vec;
    ProbabilisticAssemblingModel vj_single_model(BENCH_DATA_FOLDER + "../../models/hTRA", EMPTY);
    YMIR_BENCHMARK("VJ meta", vj_single_model.buildGraphs(cloneset_vj, SAVE_METADATA, NO_ERRORS, false))
    YMIR_BENCHMARK("VJ prob", vec = vj_single_model.computeFullProbabilities(cloneset_vj, NO_ERRORS, SUM_PROBABILITY, false))
    std::cout << loglikelihood(vec) << std::endl;

    //
    // VDJ MAAG
    //
    ProbabilisticAssemblingModel vdj_single_model(BENCH_DATA_FOLDER + "../../models/hTRB", EMPTY);
    YMIR_BENCHMARK("VDJ meta", vdj_single_model.buildGraphs(cloneset_vdj, SAVE_METADATA, NO_ERRORS, false))
    YMIR_BENCHMARK("VDJ prob", vdj_single_model.computeFullProbabilities(cloneset_vdj, NO_ERRORS, SUM_PROBABILITY, false))

    //
    // VJ inference
    //
    parser.openAndParse(BENCH_DATA_FOLDER + input_alpha_file_nonc,
                        &cloneset_vj_noncoding,
                        vj_single_genes,
                        VJ_RECOMB,
                        AlignmentColumnOptions(AlignmentColumnOptions::OVERWRITE,
                                               AlignmentColumnOptions::OVERWRITE),
                        vdj_aligner_parameters_nuc);

    YMIR_BENCHMARK("VJ EM",
                   logLvec = EMAlgorithm().statisticalInference(cloneset_vj, vj_single_model,
                                                                EMAlgorithm::AlgorithmParameters()
                                                                        .set("niter", 30),
                                                                NO_ERRORS))
//
//    YMIR_BENCHMARK("VJ SG",
//                   logLvec = SGAlgorithm().statisticalInference(cloneset_vj, vj_single_model,
//                                                                SGAlgorithm::AlgorithmParameters()
//                                                                        .set("niter", 10)
//                                                                        .set("block.size", 5000)
//                                                                        .set("alpha", .7)
//                                                                        .set("beta", 1.)
//                                                                        .set("K", 2.)
//                                                                        .set("prebuild", false)
//                                                                        .set("sample", 100000),
//                                                                NO_ERRORS))


    //
    // VDJ inference
    //
    parser.openAndParse(BENCH_DATA_FOLDER + input_beta_file_nonc,
                        &cloneset_vdj_noncoding,
                        vdj_single_genes,
                        VDJ_RECOMB,
                        AlignmentColumnOptions(AlignmentColumnOptions::OVERWRITE,
                                               AlignmentColumnOptions::OVERWRITE,
                                               AlignmentColumnOptions::OVERWRITE),
                        vdj_aligner_parameters_nuc);

    YMIR_BENCHMARK("VDJ EM",
                   logLvec = EMAlgorithm().statisticalInference(cloneset_vdj, vdj_single_model,
                                                      EMAlgorithm::AlgorithmParameters()
                                                              .set("niter", 30),
                                                                NO_ERRORS))
//
//    YMIR_BENCHMARK("VDJ SG",
//                   logLvec = SGAlgorithm().statisticalInference(cloneset_vdj, vdj_single_model,
//                                                      SGAlgorithm::AlgorithmParameters()
//                                                              .set("niter", 10)
//                                                              .set("block.size", 5000)
//                                                              .set("alpha", .7)
//                                                              .set("beta", 1.)
//                                                              .set("K", 2.)
//                                                              .set("prebuild", false)
//                                                              .set("sample", 100000),
//                                                                NO_ERRORS))

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
