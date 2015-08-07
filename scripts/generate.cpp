//
// Created by Vadim N. on 04/04/2015.
//


#include <ostream>

#include "writer.h"
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
    std::string model_path(argv[1]), out_file_path(argv[3]);
    size_t count = std::stoi(argv[2]);

    ProbabilisticAssemblingModel model(model_path);

    Cloneset gen_rep = model.generateSequences(count);

    RepertoireWriter writer;
    writer.write(out_file_path, gen_rep);

    return 0;
}