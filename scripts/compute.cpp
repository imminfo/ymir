//
// Created by Vadim N. on 04/04/2015.
//
//#include <omp.h>
//#include <stdio.h>
//
//int main() {
//#pragma omp parallel
//    printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
//}


#include <ostream>

#include "Model"


using namespace ymir;

/**
 * \brief Main function of a script for computing generation probabitilies. It has
 * a strict order of input arguments in console:
 * argv[0] - name of the script (default);
 * argv[1] - path to an input file;
 * argv[2] - path to a model;
 * argv[3] - path to an output file;
 * argv[4] - 0 if model should use stored gene usage; 1 if model should recompute gene usage from the input file.
 * argv[5] - 'n' for nucleotide sequences, 'a' for amino acid sequences.
 */
int main(int argc, char* argv[]) {
    std::string in_file_path(argv[1]),
            model_path(argv[2]),
            out_file_path(argv[3]);
    bool recompute_genes = std::stoi(argv[4]);

    std::cout << "Input file:\t" << in_file_path << std::endl;
    std::cout << "Model path:\t" << model_path << std::endl;
    std::cout << "Output file:\t" << out_file_path << std::endl;
    std::cout << std::endl;

    ProbabilisticAssemblingModel model(model_path);
    std::cout << std::endl;

    if (model.status()) {

        //
        // Nucleotide
        //
//        if (argv[5] == "n") {
        ParserNuc parser(new NaiveCDR3NucleotideAligner(model.gene_segments(), VDJAlignerParameters(3)));
        ClonesetNuc cloneset;

        if (parser.openAndParse(in_file_path,
                                &cloneset,
                                model.gene_segments(),
                                model.recombination(),
                                AlignmentColumnOptions(AlignmentColumnOptions::USE_PROVIDED,
                                                       AlignmentColumnOptions::USE_PROVIDED,
                                                       AlignmentColumnOptions::OVERWRITE))) {
            if (recompute_genes) {
                std::cout << std::endl;
                std::cout << "Recomputing gene usage on " << (size_t) cloneset.noncoding().size() << " clonotypes." << std::endl;
                model.updateGeneUsage(cloneset);
            }
            std::cout << std::endl;

            ClonesetAA cloneset_aa;
            std::vector<prob_t> prob_vec;

            if (argv[5][0] == 'n') {
                prob_vec = model.computeFullProbabilities(cloneset, NO_ERRORS);
            } else {
                std::cout << "Converting nucleotide clonotypes to amino acid clonotypes..." << std::endl;
                CDR3AminoAcidAligner(model.gene_segments(), VDJAlignerParameters(3)).toAminoAcid(cloneset, &cloneset_aa);
                std::cout << "Done." << std::endl << std::endl;

                prob_vec = model.computeFullProbabilities(cloneset_aa);
            }


            std::ofstream ofs;
            ofs.open(out_file_path);

            std::cout << std::endl;
            std::cout << "Generation probabilities statistics:" << std::endl;
            prob_summary(prob_vec);

            if (ofs.is_open()) {
                for (auto i = 0; i < prob_vec.size(); ++i) {
                    ofs << prob_vec[i] << std::endl;
                }
            } else {
                std::cout << "Problems with the output stream. Terminating..." << std::endl;
            }
            ofs.close();
        } else {
            std::cout << "Problems in parsing the input file. Terminating..." << std::endl;
        }
//        }

        //
        // Amino acid
        //
//        else {
//            ParserAA parser;
//            ClonesetAA cloneset;
//
//            if (parser.openAndParse(in_file_path,
//                                    &cloneset,
//                                    model.gene_segments(),
//                                    model.recombination(),
//                                    AlignmentColumnOptions(AlignmentColumnOptions::USE_PROVIDED,
//                                                           AlignmentColumnOptions::USE_PROVIDED,
//                                                           AlignmentColumnOptions::OVERWRITE))) {
////                if (recompute_genes) {
////                    std::cout << std::endl;
////                    std::cout << "Recomputing gene usage on " << (size_t) cloneset.noncoding().size() << " clonotypes." << std::endl;
////                    model.updateGeneUsage(cloneset);
////                }
//
//                std::cout << std::endl;
//                auto prob_vec = model.computeFullProbabilities(cloneset);
//
//                std::ofstream ofs;
//                ofs.open(out_file_path);
//
//                std::cout << std::endl;
//                std::cout << "Generation probabilities statistics:" << std::endl;
//                prob_summary(prob_vec);
//
//                if (ofs.is_open()) {
//                    for (auto i = 0; i < prob_vec.size(); ++i) {
//                        ofs << prob_vec[i] << std::endl;
//                    }
//                } else {
//                    std::cout << "Problems with the output stream. Terminating..." << std::endl;
//                }
//                ofs.close();
//            } else {
//                std::cout << "Problems in parsing the input file. Terminating..." << std::endl;
//            }
//        }

    } else {
        std::cout << "Problems with the model. Terminating..." << std::endl;
    }

    return 0;
}