//
// Created by Vadim N. on 16/03/2016.
//

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

    std::cout << "Data folder:\t[" << BENCH_DATA_FOLDER << "]" << endl;

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


    VDJRecombinationGenes vj_single_genes2 = vj_single_genes;
    VDJRecombinationGenes vj_single_genes3(vj_single_genes);
    vj_single_genes2 = vj_single_genes3;
    vj_single_genes3 = vj_single_genes;

//
    Cloneset cloneset_vj;
    YMIR_BENCHMARK("Parsing VJ",
                   parser.openAndParse(BENCH_DATA_FOLDER + input_alpha_file,
                                       &cloneset_vj,
                                       vj_single_genes,
                                       NUCLEOTIDE,
                                       VJ_RECOMB,
                                       AlignmentColumnOptions(AlignmentColumnOptions::OVERWRITE, AlignmentColumnOptions::OVERWRITE),
                                       VDJAlignerParameters(2)))

    //
    // TCR beta chain repertoire - VDJ recombination
    //
//    VDJRecombinationGenes vdj_single_genes("Vgene",
//                                    BENCH_DATA_FOLDER + "trbv.txt",
//                                    "Jgene",
//                                    BENCH_DATA_FOLDER + "trbj.txt",
//                                    "Dgene",
//                                    BENCH_DATA_FOLDER + "trbd.txt");
//
//    Cloneset cloneset_vdj;
//    YMIR_BENCHMARK("Parsing VDJ",
//                   parser.openAndParse(BENCH_DATA_FOLDER + input_beta_file,
//                                       &cloneset_vdj,
//                                       vdj_single_genes,
//                                       NUCLEOTIDE,
//                                       VDJ_RECOMB,
//                                       AlignmentColumnOptions()
//                                               .setV(AlignmentColumnOptions::USE_PROVIDED)
//                                               .setD(AlignmentColumnOptions::OVERWRITE)
//                                               .setJ(AlignmentColumnOptions::USE_PROVIDED),
//                                       VDJAlignerParameters(3)))

    //
    // VJ MAAG
    //
    ProbabilisticAssemblingModel vj_single_model(BENCH_DATA_FOLDER + "../../models/hTRA", EMPTY);

    //
    // VDJ MAAG
    //
//    ProbabilisticAssemblingModel vdj_single_model(BENCH_DATA_FOLDER + "../../models/hTRB", EMPTY);


    //
    // Inference
    //
    vector<int> vec_sample = {10000, 25000, 50000, 100000, 150000};
    vec_sample = { 100000 };
    vector<int> vec_block = {500, 1000, 2000, 5000, 10000}; //, 2000, 5000, 10000};
//    vec_block = {2000, 5000, 10000};
    vec_block = { std::stoi(argv[2]) };
    vector<double> vec_alpha = {.5, .6, .7, .8, .9};
    vec_alpha = { std::stof(argv[3]) };
    vector<double> vec_beta =  {.1, .5, 1, 3};
    vec_beta = { std::stof(argv[4]) };
    vector<double> vec_K =     {1, 2, 3};
    vec_K = { std::stod(argv[5]) };
    ErrorMode error_mode = COMPUTE_ERRORS;

//    int niter, sample, block;
//    double alpha;
//    niter = 10;
//    int sample = 10000;

//    RUN_EM_INFERENCE(string("vj"), cloneset_vj, vj_single_model, niter, sample, error_mode)
//    RUN_EM_INFERENCE(string("vdj"), cloneset_vdj, vdj_single_model, niter, sample, error_mode)
//    for(auto val_sample: vec_sample) {
//        RUN_EM_INFERENCE(string("vj"), cloneset_vj, vj_single_model, 40, val_sample, error_mode)
//    }

//    ModelParameterVector new_param_vec = vj_single_model.event_probabilities();
//    new_param_vec.fill(1);
//    new_param_vec.normaliseEventFamilies();
//    vj_single_model.updateModelParameterVector(new_param_vec);
//
//    auto rep_nonc = cloneset_vj.noncoding().sample(2500);
//    auto maag_rep = vj_single_model.buildGraphs(rep_nonc, SAVE_METADATA, error_mode, NUCLEOTIDE, true);
//    vector<prob_t> prob_vec(maag_rep.size(), 0);
//
//    vector<bool> good_clonotypes(maag_rep.size(), true);
//    for (size_t i = 0; i < maag_rep.size(); ++i) {
//        if (rep_nonc[i].is_good()) {
//            prob_vec[i] = maag_rep[i].fullProbability();
//            if (std::isnan(prob_vec[i]) || prob_vec[i] == 0) {
//                good_clonotypes[i] = false;
//            }
//        } else {
//            good_clonotypes[i] = false;
//        }
//    }

    for(auto val_sample: vec_sample) {
        for(auto val_block: vec_block) {
            for (auto val_alpha: vec_alpha) {
                for (auto val_beta: vec_beta) {
                    for (auto val_K: vec_K) {
                        RUN_SG_INFERENCE(string("vj"), cloneset_vj, vj_single_model, 40, val_block, val_alpha, val_beta, val_K, val_sample, error_mode)
                    }
                }
            }
        }
    }

}
