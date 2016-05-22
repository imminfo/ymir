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

    time_t vj_single_parse, vdj_single_parse,
            vj_single_prob, vj_single_meta,
            vj_single_infer, vdj_single_prob,
            vdj_single_meta, vdj_single_infer;

    time_t vj_good_parse, vdj_good_parse,
            vj_good_prob, vj_good_meta,
            vj_good_infer, vdj_good_prob,
            vdj_good_meta, vdj_good_infer;
    
    std::string BENCH_DATA_FOLDER = argv[1];


    ParserNuc parser(new NaiveCDR3NucleotideAligner());

    string input_alpha_file = "alpha.100k.txt";
    string input_beta_file = "beta.100k.txt";


    //
    // TCR alpha chain repertoire - VJ recombination
    //
    VDJRecombinationGenes vj_single_genes("Vgene",
                                   BENCH_DATA_FOLDER + "trav.txt",
                                   "Jgene",
                                   BENCH_DATA_FOLDER + "traj.txt");

    ClonesetNuc cloneset_vj;
    YMIR_BENCHMARK("Parsing VJ",
                   parser.openAndParse(BENCH_DATA_FOLDER + input_alpha_file,
                                       &cloneset_vj,
                                       vj_single_genes,
                                       VJ_RECOMB,
                                       AlignmentColumnOptions(AlignmentColumnOptions::OVERWRITE,
                                                              AlignmentColumnOptions::OVERWRITE),
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

    ClonesetNuc cloneset_vdj;
    YMIR_BENCHMARK("Parsing VDJ",
                   parser.openAndParse(BENCH_DATA_FOLDER + input_beta_file,
                                       &cloneset_vdj,
                                       vdj_single_genes,
                                       VDJ_RECOMB,
                                       AlignmentColumnOptions(AlignmentColumnOptions::OVERWRITE,
                                                              AlignmentColumnOptions::OVERWRITE,
                                                              AlignmentColumnOptions::OVERWRITE),
                                       VDJAlignerParameters(3)))

    //
    // VJ MAAG
    //
    ProbabilisticAssemblingModel vj_single_model(BENCH_DATA_FOLDER + "../../models/hTRA", EMPTY);
//    YMIR_BENCHMARK("VJ meta", vj_single_model.buildGraphs(cloneset_vj, SAVE_METADATA, NO_ERRORS))
//    YMIR_BENCHMARK("VJ prob", vj_single_model.computeFullProbabilities(cloneset_vj, NO_ERRORS))


    //
    // VDJ MAAG
    //
//    ProbabilisticAssemblingModel vdj_single_model(BENCH_DATA_FOLDER + "../../models/hTRB", EMPTY);

//    YMIR_BENCHMARK("VDJ meta", vdj_single_model.buildGraphs(cloneset_vdj, SAVE_METADATA, NO_ERRORS))
//    YMIR_BENCHMARK("VDJ prob", vdj_single_model.computeFullProbabilities(cloneset_vdj, NO_ERRORS))


    //
    // VJ inference
    //
    YMIR_BENCHMARK("VJ EM",
                   logLvec = EMAlgorithm().statisticalInference(cloneset_vj, vj_single_model,
                                                                EMAlgorithm::AlgorithmParameters()
                                                                        .set("niter", 10)
                                                                        .set("sample", 100000),
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
//    YMIR_BENCHMARK("VDJ EM",
//                   logLvec = EMAlgorithm().statisticalInference(cloneset_vdj, vdj_single_model,
//                                                      EMAlgorithm::AlgorithmParameters()
//                                                              .set("niter", 10)
//                                                              .set("sample", 100000),
//                                                                NO_ERRORS))
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
