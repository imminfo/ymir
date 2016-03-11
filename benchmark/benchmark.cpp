//
// Created by Vadim N. on 18/03/2015.
//

#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_


#define YMIR_BENCHMARK(expr, time_var) { time_var = std::chrono::system_clock::now(); \ }


#include <chrono>
#include <ctime>

#include "Inference"


using namespace ymir;


int main(int argc, char* argv[]) {
    std::chrono::system_clock::time_point tp1, tp2;
    
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


    //
    // TCR alpha chain repertoire - VJ recombination
    //
    VDJRecombinationGenes vj_single_genes("Vgene",
                                   BENCH_DATA_FOLDER + "trav.txt",
                                   "Jgene",
                                   BENCH_DATA_FOLDER + "traj.txt");

    tp1 = std::chrono::system_clock::now();
    Cloneset cloneset_vj;
    parser.openAndParse(BENCH_DATA_FOLDER + "alpha.250k.txt",
                        &cloneset_vj,
                        vj_single_genes,
                        NUCLEOTIDE,
                        VJ_RECOMB,
                        AlignmentColumnOptions()
                                .setV(AlignmentColumnOptions::OVERWRITE)
                                .setJ(AlignmentColumnOptions::OVERWRITE),
                        VDJAlignerParameters(2));

    tp2 = std::chrono::system_clock::now();
    vj_single_parse = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);

    //
    // TCR beta chain repertoire - VDJ recombination
    //
    VDJRecombinationGenes vdj_single_genes("Vgene",
                                    BENCH_DATA_FOLDER + "trbv.txt",
                                    "Jgene",
                                    BENCH_DATA_FOLDER + "trbj.txt",
                                    "Dgene",
                                    BENCH_DATA_FOLDER + "trbd.txt");

    tp1 = std::chrono::system_clock::now();
    Cloneset cloneset_vdj;
//    parser.openAndParse(BENCH_DATA_FOLDER + "beta.250k.txt",
//                        &cloneset_vdj,
//                        vdj_single_genes,
//                        NUCLEOTIDE,
//                        VDJ_RECOMB,
//                        AlignmentColumnOptions()
//                                .setV(AlignmentColumnOptions::USE_PROVIDED)
//                                .setJ(AlignmentColumnOptions::USE_PROVIDED)
//                                .setD(AlignmentColumnOptions::OVERWRITE),
//                        VDJAlignerParameters(1, 3));

    tp2 = std::chrono::system_clock::now();
    vdj_single_parse = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);

    //
    // VJ MAAG
    //
    ProbabilisticAssemblingModel vj_single_model(BENCH_DATA_FOLDER + "../../models/hTRA"); //, EMPTY);

    tp1 = std::chrono::system_clock::now();
//    vj_single_model.buildGraphs(cloneset_vj, SAVE_METADATA, NO_ERRORS);
    tp2 = std::chrono::system_clock::now();
    vj_single_meta = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);

    tp1 = std::chrono::system_clock::now();
//    vj_single_model.computeFullProbabilities(cloneset_vj, NO_METADATA, NO_ERRORS);
    tp2 = std::chrono::system_clock::now();
    vj_single_prob = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);


    //
    // VDJ MAAG
    //
    ProbabilisticAssemblingModel vdj_single_model(BENCH_DATA_FOLDER + "../../models/hTRB", EMPTY);

    tp1 = std::chrono::system_clock::now();
//    vdj_single_model.buildGraphs(cloneset_vdj, SAVE_METADATA, NO_ERRORS);
    tp2 = std::chrono::system_clock::now();
    vdj_single_meta = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);

    tp1 = std::chrono::system_clock::now();
//    vdj_single_model.computeFullProbabilities(cloneset_vdj, NO_METADATA, NO_ERRORS);
    tp2 = std::chrono::system_clock::now();
    vdj_single_prob = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);


    //
    // VJ inference
    //
    tp1 = std::chrono::system_clock::now();
    EMAlgorithm().statisticalInference(cloneset_vj,
                                       vj_single_model,
                                       EMAlgorithm::AlgorithmParameters().set("niter", 50),
                                       COMPUTE_ERRORS);
//                                       COMPUTE_ERRORS);
//    BlockEMAlgorithm().statisticalInference(cloneset_vj,
//                                       vj_single_model,
//                                       EMAlgorithm::AlgorithmParameters().set("niter", 30).set("block.size", 2000),
//                                       NO_ERRORS);
    tp2 = std::chrono::system_clock::now();
    vj_single_infer = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);


    //
    // VDJ inference
    //
    tp1 = std::chrono::system_clock::now();
//    EMAlgorithm().statisticalInference(cloneset_vdj,
//                                       vdj_single_model,
//                                       EMAlgorithm::AlgorithmParameters().set("niter", 10));
    tp2 = std::chrono::system_clock::now();
    vdj_single_infer = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);


    //
    // Results
    //
    cout << "========================" << endl << "Results:" << endl;

    cout << endl;
    cout << "TCR alpha repertoire:\tsingle V-J alignment" << endl;
    cout << "TCR beta repertoire:\tsingle V-J alignment, all D alignments" << endl;

    cout << "Parsing VJ, seconds:\t" << vj_single_parse << endl;
    cout << "Parsing VDJ, seconds:\t" << vdj_single_parse << endl;

    cout << "VJ MAAG computing, seconds:\t" << vj_single_prob << endl;
    cout << "VJ MAAG metadata, seconds:\t" << vj_single_meta << endl;
    cout << "VJ inference 10 iter, seconds:\t" << vj_single_infer << endl;

    cout << "VDJ MAAG computing, seconds:\t" << vdj_single_prob << endl;
    cout << "VDJ MAAG metadata, seconds:\t" << vdj_single_meta << endl;
    cout << "VDJ inference 10 iter, seconds:\t" << vdj_single_infer << endl;

    cout << endl;
    cout << "TCR alpha repertoire:\tall good V-J alignments" << endl;
    cout << "TCR beta repertoire:\tall good V-J alignments, all D alignments" << endl;

    cout << "Parsing VJ, seconds:\t" << vj_good_parse << endl;
    cout << "Parsing VDJ, seconds:\t" << vdj_good_parse << endl;

    cout << "VJ MAAG computing, seconds:\t" << vj_good_prob << endl;
    cout << "VJ MAAG metadata, seconds:\t" << vj_good_meta << endl;
    cout << "VJ inference 10 iter, seconds:\t" << vj_good_infer << endl;

    cout << "VDJ MAAG computing, seconds:\t" << vdj_good_prob << endl;
    cout << "VDJ MAAG metadata, seconds:\t" << vdj_good_meta << endl;
    cout << "VDJ inference 10 iter, seconds:\t" << vdj_good_infer << endl;

    return 0;
}

#endif //_BENCHMARK_H_
