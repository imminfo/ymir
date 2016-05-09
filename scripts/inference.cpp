//
// Created by Vadim N. on 04/04/2015.
//


#include <ostream>

#include "Inference"


using namespace ymir;

/**
 * \brief Main function of a script for the statistical inference of marginal probabilities. It has
 * a strict order of input arguments in console:
 * argv[0] - name of the script (default);
 * argv[1] - path to a file with immune receptors
 * argv[2] - path to a basic model
 * argv[3] - path to an output model folder
 * argv[4] - which algorithm to use
 * argv[5], argv[6], ... - algorithms parameters in a form "name" - "value"
 */
int main(int argc, char* argv[]) {
    std::string in_file_path(argv[1]),
            model_path(argv[2]),
            out_model_path(argv[3]),
            algo_id(argv[4]), algo_name = "noname";

    StatisticalInferenceAlgorithm *infer_algo = nullptr;
    StatisticalInferenceAlgorithm::AlgorithmParameters params;
    if (algo_id == "em") {
        infer_algo = new EMAlgorithm();
        algo_name = "EM-algorithm";
        params = params.set("niter", std::stoi(argv[6]));  // argv[5] == "niter"
    }

    std::cout << "Input cloneset file:\t" << in_file_path << std::endl;
    std::cout << "Basic model path:\t" << model_path << std::endl;
    std::cout << "Output model path:\t" << out_model_path << std::endl;
    std::cout << "Algorithm:\t" << algo_name << std::endl;
    std::cout << std::endl;

    ProbabilisticAssemblingModel model(model_path);
    std::cout << std::endl;

    if (model.status()) {
        ParserNuc parser;
        ClonesetNuc cloneset;

        if (parser.openAndParse( in_file_path,
                                 &cloneset,
                                 model.gene_segments(),
                                 model.recombination(),
                                 AlignmentColumnOptions(AlignmentColumnOptions::OVERWRITE,
                                                        AlignmentColumnOptions::OVERWRITE,
                                                        AlignmentColumnOptions::OVERWRITE))) {

            std::cout << std::endl;
            infer_algo->statisticalInference(cloneset, model, params);

            if (model.save(out_model_path)) {

            } else {
                std::cout << "Problems in saving the resulting model. Terminating..." << std::endl;
            }
        } else {
            std::cout << "Problems in parsing the input file. Terminating..." << std::endl;
        }
    } else {
        std::cout << "Problems with the model. Terminating..." << std::endl;
    }

    if (infer_algo) {
        delete infer_algo;
    }
    return 0;
}