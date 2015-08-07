//
// Created by Vadim N. on 04/04/2015.
//


#include <ostream>

#include "probabilisticassemblingmodel.h"


using namespace ymir;

/**
 * \brief Main function of a script for generation of pre-selected receptors. It has
 * a strict order of input arguments in console:
 * argv[0] - name of the script (default);
 * argv[1] - path to a model
 * argv[2] - number of immune receptors to generate
 * argv[3] - path to an output file
 */
int main(int argc, char* argv[]) {
    std::string in_file_path(argv[1]), model_path(argv[2]), out_file_path(argv[3]);
    bool recompute_genes = std::stoi(argv[4]);

    ProbabilisticAssemblingModel model(model_path);

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
    for (auto i = 0; i < prob_vec.size(); ++i) {
        ofs << prob_vec[i] << std::endl;
    }
    ofs.close();

    return 0;
}