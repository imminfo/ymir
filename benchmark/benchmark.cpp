//
// Created by Vadim N. on 18/03/2015.
//

#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#define BENCH_DATA_FOLDER string("/Users/vdn/ymir/benchmark/data/")


#include <chrono>
#include <ctime>

#include "parser.h"
#include "statisticalinferencealgorithm.h"

using namespace ymir;


int main() {
    std::chrono::system_clock::time_point tp1, tp2;
    time_t vj_prob, vj_meta, vdj_prob, vdj_meta;

    RepertoireParser parser;
    parser.loadConfig(BENCH_DATA_FOLDER + "../../parsers/mitcrdots.json");


    //
    // TCR alpha chain repertoire - VJ recombination
    //
    VDJRecombinationGenes vj_genes("Vgene",
                                   BENCH_DATA_FOLDER + "trav.txt",
                                   "Jgene",
                                   BENCH_DATA_FOLDER + "traj.txt");

    tp1 = std::chrono::system_clock::now();
    Cloneset cloneset_vj;
    parser.parse(BENCH_DATA_FOLDER + "mitcr.alpha.500k.txt",
                 &cloneset_vj,
                 vj_genes,
                 RepertoireParser::AlignmentColumnOptions()
                         .setV(RepertoireParser::MAKE_IF_NOT_FOUND)
                         .setJ(RepertoireParser::MAKE_IF_NOT_FOUND)
                         .setD(RepertoireParser::SKIP));
    tp2 = std::chrono::system_clock::now();
    cout << "Parsing VJ, seconds: " << (std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1)) << endl;

    //
    // TCR beta chain repertoire - VDJ recombination
    //
    VDJRecombinationGenes vdj_genes("Vgene",
                                    BENCH_DATA_FOLDER + "trbv.txt",
                                    "Jgene",
                                    BENCH_DATA_FOLDER + "trbj.txt",
                                    "Dgene",
                                    BENCH_DATA_FOLDER + "trbd.txt");

    tp1 = std::chrono::system_clock::now();
    Cloneset cloneset_vdj;
    parser.parse(BENCH_DATA_FOLDER + "mitcr.beta.500k.txt",
                 &cloneset_vdj,
                 vdj_genes,
                 RepertoireParser::AlignmentColumnOptions()
                         .setV(RepertoireParser::MAKE_IF_NOT_FOUND)
                         .setJ(RepertoireParser::MAKE_IF_NOT_FOUND)
                         .setD(RepertoireParser::OVERWRITE));
    tp2 = std::chrono::system_clock::now();
    cout << "Parsing VDJ, seconds: " << (std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1)) << endl;

    //
    // VJ MAAG
    //
    ProbabilisticAssemblingModel vj_model(BENCH_DATA_FOLDER + "../../models/hTRA");

    tp1 = std::chrono::system_clock::now();
    vj_model.buildGraphs(cloneset_vj, true);
    tp2 = std::chrono::system_clock::now();
    vj_meta = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);

    tp1 = std::chrono::system_clock::now();
    vj_model.computeFullProbabilities(cloneset_vj);
    tp2 = std::chrono::system_clock::now();
    vj_prob = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);


    //
    // VDJ MAAG
    //
    ProbabilisticAssemblingModel vdj_model(BENCH_DATA_FOLDER + "../../models/hTRB");

    tp1 = std::chrono::system_clock::now();
    MAAGRepertoire(vdj_model.buildGraphs(cloneset_vdj, true));
    tp2 = std::chrono::system_clock::now();
    vdj_meta = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);

    tp1 = std::chrono::system_clock::now();
    vdj_model.computeFullProbabilities(cloneset_vdj);
    tp2 = std::chrono::system_clock::now();
    vdj_prob = std::chrono::system_clock::to_time_t(tp2)- std::chrono::system_clock::to_time_t(tp1);

    //
    // Results
    //
    cout << "========================" << endl << "Results:" << endl;
    cout << "VJ MAAG computing, seconds: " << vj_prob << endl;
    cout << "VJ MAAG metadata, seconds: " << vj_meta << endl;
    cout << "VDJ MAAG computing, seconds: " << vdj_prob << endl;
    cout << "VDJ MAAG metadata, seconds: " << vdj_meta << endl;


//    Cloneset repertoire50K, repertoire200K, repertoire500K;
//  With and without parallelisation

    // Compute simple one-alignment full probabilities.
    // 1 V ; 1 J
    //    repertoire1K;
//    repertoire10K;
//    repertoire100K;
    // 1000
    // 10000
    // 500000

    // 1 V ; 2 Ds ; 1 J
    //    repertoire1K;
//    repertoire10K;
//    repertoire100K;
    // 1000
    // 10000
    // 500000


    // Compute simple one-alignment full probabilities with storing additional data like CDR3 sequences.
    // 1 V ; 1 J
    // 1000
    // 10000
    // 500000

    // 1 V ; 2 Ds ; 1 J
    // 1000
    // 10000
    // 500000


    // Compute various multi-alignment full probabilities.
    // 2 Vs ; 1 Js
    // 1000
    // 10000
    // 500000

    // 5 Vs ; 1 Js
    // 1000
    // 10000
    // 500000

    // 5 Vs ; 5 Js
    // 1000
    // 10000
    // 500000

    // 2 Vs ; 2 Ds ; 1 Js
    // 1000
    // 10000
    // 500000

    // 5 Vs ; 2 Ds ; 1 Js
    // 1000
    // 10000
    // 500000

    // 5 Vs ; 2 Ds ; 5 Js
    // 1000
    // 10000
    // 500000


    // Compute various multi-alignment full probabilities with full graph building.
    // 2 Vs ; 1 Js
    // 1000
    // 10000
    // 500000

    // 5 Vs ; 1 Js
    // 1000
    // 10000
    // 500000

    // 5 Vs ; 5 Js
    // 1000
    // 10000
    // 500000

    // 2 Vs ; 2 Ds ; 1 Js
    // 1000
    // 10000
    // 500000

    // 5 Vs ; 2 Ds ; 1 Js
    // 1000
    // 10000
    // 500000

    // 5 Vs ; 2 Ds ; 5 Js
    // 1000
    // 10000
    // 500000



    return 0;
}

#endif //_BENCHMARK_H_
