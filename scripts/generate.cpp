//
// Created by Vadim N. on 04/04/2015.
//


#include <ostream>

#include "Model"


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

    std::cout << "Model path:\t" << model_path << std::endl;
    std::cout << "Output file:\t" << out_file_path << std::endl;
    std::cout << "Number of clonotypes:\t" << (size_t) count << std::endl;
    std::cout << std::endl;

    ProbabilisticAssemblingModel model(model_path);

    if (model.status()) {
        std::cout << std::endl;

        size_t to_generate, generated = 0;
        size_t block_size = 100000;
        RepertoireWriter writer;
        while (generated < count) {
            to_generate = std::min(count - generated, block_size);
            generated += to_generate;
            Cloneset gen_rep = model.generateSequences(to_generate, false);
            if (!writer.write(out_file_path, gen_rep, model.gene_segments(), true)) {
                std::cout << "Problems in writing the output file. Terminating..." << std::endl;
            }
            std::cout << "Generated " << (size_t) generated << "/" << (size_t) count << std::endl;
        }
    } else {
        std::cout << "Problems with the model. Terminating..." << std::endl;
    }

    return 0;
}