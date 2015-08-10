//
// Created by Vadim N. on 07/08/2015.
//

#ifndef YMIR_WRITER_H
#define YMIR_WRITER_H


#include <ostream>

#include "repertoire.h"
#include "genesegment.h"


namespace ymir {

    class RepertoireWriter {
    public:


        RepertoireWriter() { }


        virtual ~RepertoireWriter() { }


        bool write(const std::string &filepath,
                   const ClonesetView &cloneset,
                   const VDJRecombinationGenes& gene_segments) const {
            std::ofstream ofs;

            ofs.open(filepath);

            if (ofs.is_open()) {

                // write the header
                ofs << "Nucleotide sequence" << '\t' <<
                        "Amino acid sequence" << '\t' <<
                        "Variable" << '\t' <<
                        "Diversity" << '\t' <<
                        "Joining" << '\t' <<
                        "V end" << '\t' <<
                        "D start" << '\t' <<
                        "D end" << '\t' <<
                        "J start" << std::endl;

                // write clonotypes
                for (auto i = 0; i < cloneset.size(); ++i) {
                    if (cloneset[i].is_nucleotide()) {
                        ofs << cloneset[i].sequence() << '\t';
                        ofs << translate(cloneset[i].sequence()) << '\t';
                    } else {
                        ofs << "" << '\t';
                        ofs << cloneset[i].sequence() << '\t';
                    }

                    for (auto seg_i = 0; seg_i < cloneset[i].nVar(); ++seg_i) {
                        ofs << gene_segments.V()[cloneset[i].getVar(seg_i)].allele;
                        if (seg_i < cloneset[i].nVar() - 1) {
                            ofs << ',';
                        } else {
                            ofs << '\t';
                        }
                    }

                    ofs << "\t\t";

                    for (auto seg_i = 0; seg_i < cloneset[i].nJoi(); ++seg_i) {
                        ofs << gene_segments.J()[cloneset[i].getJoi(seg_i)].allele;
                        if (seg_i < cloneset[i].nJoi() - 1) {
                            ofs << ',';
                        } else {
                            ofs << '\t';
                        }
                    }
                }

                ofs.close();
                return true;
            }

            return false;
        }

    protected:


    };
}


#endif //YMIR_WRITER_H
