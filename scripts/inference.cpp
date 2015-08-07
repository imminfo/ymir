//
// Created by Vadim N. on 04/04/2015.
//


#include <ostream>

#include "statisticalinferencealgorithm.h"


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
            algo_id(argv[4]);

    ProbabilisticAssemblingModel model(model_path);
    StatisticalInferenceAlgorithm *infer_algo;
    StatisticalInferenceAlgorithm::AlgorithmParameters params;

    if (algo_id == "em") {
        params = params.set("niter", std::stoi(argv[6]));  // argv[5] == "niter"
    }

    RepertoireParser parser;
    Cloneset cloneset;
    parser.parse(in_file_path, &cloneset, model.gene_segments(),
                 RepertoireParser::AlignmentColumnOptions()
                         .setV(RepertoireParser::MAKE_IF_NOT_FOUND)
                         .setJ(RepertoireParser::MAKE_IF_NOT_FOUND)
                         .setD(RepertoireParser::OVERWRITE));

    infer_algo->statisticalInference(cloneset, model, params);

    model.save(out_model_path);

    return 0;
}