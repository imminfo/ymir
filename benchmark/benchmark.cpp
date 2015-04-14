//
// Created by Vadim N. on 18/03/2015.
//

#ifndef _BENCHMARK_H_
#define _BENCHMARK_H_

#define BENCH_DATA_FOLDER string("/Users/vdn/ymir/benchmark/data/")

#include "parser.h"
#include "statisticalinferencealgorithm.h"

using namespace ymir;


int main() {

//    RepertoireParser parser;
//    parser.loadConfig(BENCH_DATA_FOLDER + "../../parsers/mitcrdots.json");
//
//    VDJRecombinationGenes genes("Vgene",
//                                BENCH_DATA_FOLDER + "trav.txt",
//                                "Jgene",
//                                BENCH_DATA_FOLDER + "traj.txt");
//
//    Cloneset cr;
//    parser.parse(BENCH_DATA_FOLDER + "mitcr.alpha.300k.txt",
//                 &cr,
//                 genes,
//                 RepertoireParser::AlignmentColumnOptions()
//                         .setV(RepertoireParser::MAKE_IF_NOT_FOUND)
//                         .setJ(RepertoireParser::MAKE_IF_NOT_FOUND)
//                         .setD(RepertoireParser::SKIP));

    RepertoireParser parser;
    parser.loadConfig(BENCH_DATA_FOLDER + "../../parsers/mitcrdots.json");

    VDJRecombinationGenes genes("Vgene",
                                BENCH_DATA_FOLDER + "trbv.txt",
                                "Jgene",
                                BENCH_DATA_FOLDER + "trbj.txt",
                                "Dgene",
                                BENCH_DATA_FOLDER + "trbd.txt");

    Cloneset cr;
    parser.parse(BENCH_DATA_FOLDER + "mitcr.beta.200k.txt",
                 &cr,
                 genes,
                 RepertoireParser::AlignmentColumnOptions()
                         .setV(RepertoireParser::MAKE_IF_NOT_FOUND)
                         .setJ(RepertoireParser::MAKE_IF_NOT_FOUND)
                         .setD(RepertoireParser::OVERWRITE));

//    Cloneset repertoire1K, repertoire10K, repertoire100K;
//  With and without parallelisation

    // Compute simple one-alignment full probabilities.
    // 1 V ; 1 J
    //    repertoire1K;
//    repertoire10K;
//    repertoire100K;
    // 1000
    // 10000
    // 100000

    // 1 V ; 2 Ds ; 1 J
    //    repertoire1K;
//    repertoire10K;
//    repertoire100K;
    // 1000
    // 10000
    // 100000


    // Compute simple one-alignment full probabilities with storing additional data like CDR3 sequences.
    // 1 V ; 1 J
    // 1000
    // 10000
    // 100000

    // 1 V ; 2 Ds ; 1 J
    // 1000
    // 10000
    // 100000


    // Compute various multi-alignment full probabilities.
    // 2 Vs ; 1 Js
    // 1000
    // 10000
    // 100000

    // 5 Vs ; 1 Js
    // 1000
    // 10000
    // 100000

    // 5 Vs ; 5 Js
    // 1000
    // 10000
    // 100000

    // 2 Vs ; 2 Ds ; 1 Js
    // 1000
    // 10000
    // 100000

    // 5 Vs ; 2 Ds ; 1 Js
    // 1000
    // 10000
    // 100000

    // 5 Vs ; 2 Ds ; 5 Js
    // 1000
    // 10000
    // 100000


    // Compute various multi-alignment full probabilities with full graph building.
    // 2 Vs ; 1 Js
    // 1000
    // 10000
    // 100000

    // 5 Vs ; 1 Js
    // 1000
    // 10000
    // 100000

    // 5 Vs ; 5 Js
    // 1000
    // 10000
    // 100000

    // 2 Vs ; 2 Ds ; 1 Js
    // 1000
    // 10000
    // 100000

    // 5 Vs ; 2 Ds ; 1 Js
    // 1000
    // 10000
    // 100000

    // 5 Vs ; 2 Ds ; 5 Js
    // 1000
    // 10000
    // 100000



    return 0;
}

#endif //_BENCHMARK_H_
