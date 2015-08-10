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

#include "parser.h"
#include "probabilisticassemblingmodel.h"


using namespace ymir;

/**
 * \brief Main function of a script for computing generation probabitilies. It has
 * a strict order of input arguments in console:
 * argv[0] - name of the script (default);
 * argv[1] - path to an input file;
 * argv[2] - path to a model;
 * argv[3] - path to an output file;
 * argv[4] - 0 if model should use stored gene usage; 1 if model should recompute gene usage from the input file.
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

    if (model.status()) {
        RepertoireParser parser;
        Cloneset cloneset;
        parser.parse(in_file_path, &cloneset, model.gene_segments(),
                     RepertoireParser::AlignmentColumnOptions()
                             .setV(RepertoireParser::MAKE_IF_NOT_FOUND)
                             .setJ(RepertoireParser::MAKE_IF_NOT_FOUND)
                             .setD(RepertoireParser::OVERWRITE));

        if (recompute_genes) {
            model.updateGeneUsage(cloneset);
        }
        auto prob_vec = model.computeFullProbabilities(cloneset);

        std::ofstream ofs;
        ofs.open(out_file_path);
        if (ofs.is_open()) {
            for (auto i = 0; i < prob_vec.size(); ++i) {
                ofs << prob_vec[i] << std::endl;
            }
            ofs.close();
        } else {
            std::cout << "Problems with the output stream. Terminating..." << std::endl;
        }
    } else {
        std::cout << "Problems with the model. Terminating..." << std::endl;
    }

    return 0;
}