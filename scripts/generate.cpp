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
 * argv[1] - path to a model;
 * argv[2] - number of immune receptors to generate;
 * argv[3] - path to an output file;
 * argv[4] - output sequence type: 0 for all, 1 for coding, 2 for noncoding.
 */
int main(int argc, char* argv[]) {
    std::string model_path(argv[1]), out_file_path(argv[3]);
    size_t count = std::stoi(argv[2]);
    SequenceCodingType coding_type = SequenceCodingType::ALL;
    if (argc >= 5) {
        switch(std::stoi(argv[4])) {
            case 0:
                std::cout << "Generate all sequence types" << std::endl;
                coding_type = SequenceCodingType::ALL;
                break;
            case 1:
                std::cout << "Generate only coding sequences" << std::endl;
                coding_type = SequenceCodingType::CODING;
                break;
            case 2:
                std::cout << "Generate only noncoding sequences" << std::endl;
                coding_type = SequenceCodingType::NONCODING;
                break;
            case 3:
                std::cout << "Generate only out-of-frame sequences" << std::endl;
                coding_type = SequenceCodingType::OUTOFFRAME;
                break;
            case 4:
                std::cout << "Generate only sequences with stop codon" << std::endl;
                coding_type = SequenceCodingType::STOPCODON;
                break;
            default:
                std::cout << "Default: generate all sequence types" << std::endl;
                coding_type = SequenceCodingType::ALL;
        }
    }

    std::cout << "Model path:\t" << model_path << std::endl;
    std::cout << "Output file:\t" << out_file_path << std::endl;
    std::cout << "Number of clonotypes:\t" << (size_t) count << std::endl;
    std::cout << std::endl;

    ProbabilisticAssemblingModel model(model_path);

    if (model.status()) {
        std::cout << std::endl;

        size_t to_generate, generated = 0, generated_coding = 0, generated_noncoding = 0;
        size_t block_size = 100000;
        RepertoireWriter writer;
        bool first_iter = false;
        while (generated < count) {
            to_generate = std::min(count - generated, block_size);
            ClonesetViewNuc gen_rep = model.generateSequences(std::max((size_t) 5000, to_generate));
            if (coding_type == SequenceCodingType::ALL) {
                if (gen_rep.size() > to_generate) { gen_rep = gen_rep.head(to_generate); }
                generated_coding +=    gen_rep.coding().size();
                generated_noncoding += gen_rep.noncoding().size();
            } else if (coding_type == SequenceCodingType::CODING) {
                gen_rep = gen_rep.coding();
                if (gen_rep.size() > to_generate) { gen_rep = gen_rep.head(to_generate); }
                generated_coding += gen_rep.size();
            } else if (coding_type == SequenceCodingType::NONCODING) {
                gen_rep = gen_rep.noncoding();
                if (gen_rep.size() > to_generate) { gen_rep = gen_rep.head(to_generate); }
                generated_noncoding += gen_rep.size();
            } else if (coding_type == SequenceCodingType::OUTOFFRAME) {
                gen_rep = gen_rep.outofframes();
                if (gen_rep.size() > to_generate) { gen_rep = gen_rep.head(to_generate); }
                generated_noncoding += gen_rep.size();
            } else if (coding_type == SequenceCodingType::STOPCODON) {
                gen_rep = gen_rep.withstopcodons();
                if (gen_rep.size() > to_generate) { gen_rep = gen_rep.head(to_generate); }
                generated_noncoding += gen_rep.size();
            } else {
                std::cout << "Something strange happened: no such sequence type." << std::endl;
            }

            generated = generated_coding + generated_noncoding;
            if (!writer.write(out_file_path, gen_rep, model.gene_segments(), first_iter)) {
                std::cout << "Problems occurred while writing the output file. Terminating..." << std::endl;
                break;
            }
            std::cout << "Generated " << (size_t) generated << "/" << (size_t) count << std::endl;
            std::cout << "\t- " << (size_t) generated_coding << " coding sequences" << std::endl;
            std::cout << "\t- " << (size_t) generated_noncoding << " noncoding sequences" << std::endl;

            first_iter = true;
        }
    } else {
        std::cout << "Problems with the model. Terminating..." << std::endl;
    }

    return 0;
}