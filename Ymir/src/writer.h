//
// Created by Vadim N. on 07/08/2015.
//

#ifndef YMIR_WRITER_H
#define YMIR_WRITER_H


#include <ostream>

#include "repertoire.h"
#include "genesegment.h"


namespace ymir {

    #define CELL_FILL(i, limit, ofs, if_letter, else_letter) { if (i < limit - 1) { ofs << if_letter; } else { ofs << else_letter; } }

    class RepertoireWriter {
    public:


        RepertoireWriter() { }


        virtual ~RepertoireWriter() { }


        template <typename ClonotypeType>
        bool write(const std::string &filepath,
                   const ClonesetView<ClonotypeType> &cloneset,
                   const VDJRecombinationGenes& gene_segments,
                   bool append = false) const {
            std::ofstream ofs;

            if (append) {
                ofs.open(filepath, std::ofstream::out | std::ofstream::app);
            } else {
                ofs.open(filepath, std::ofstream::out);
            }

            if (ofs.is_open()) {

                if (!append) {
                    // write the header
                    ofs <<  "Nucleotide sequence" << '\t' <<
                            "Amino acid sequence" << '\t' <<
                            "Variable" << '\t' <<
                            "Diversity" << '\t' <<
                            "Joining" << '\t' <<
                            "V alignment" << '\t' <<
                            "D alignment" << '\t' <<
                            "J alignment" << std::endl;
                }

                // write clonotypes
                for (auto i = 0; i < cloneset.size(); ++i) {
                    if (ClonotypeType::sequence_type == NUCLEOTIDE) {
                        ofs << cloneset[i].sequence() << '\t';
                        ofs << translate(cloneset[i].sequence()) << '\t';
                    } else {
                        ofs << "" << '\t';
                        ofs << cloneset[i].sequence() << '\t';
                    }

                    for (auto seg_i = 0; seg_i < cloneset[i].nVar(); ++seg_i) {
                        ofs << gene_segments.V()[cloneset[i].getVar(seg_i)].allele;
                        CELL_FILL(seg_i, cloneset[i].nVar(), ofs, ',', '\t')
                    }

                    if (gene_segments.is_vdj()) {
                        for (auto seg_i = 0; seg_i < cloneset[i].nDiv(); ++seg_i) {
                            ofs << gene_segments.D()[cloneset[i].getDiv(seg_i)].allele;
                            CELL_FILL(seg_i, cloneset[i].nDiv(), ofs, ',', '\t')
                        }
                    } else {
                        ofs << "\t";
                    }

                    for (auto seg_i = 0; seg_i < cloneset[i].nJoi(); ++seg_i) {
                        ofs << gene_segments.J()[cloneset[i].getJoi(seg_i)].allele;
                        CELL_FILL(seg_i, cloneset[i].nJoi(), ofs, ',', '\t')
                    }

                    for (auto seg_i = 0; seg_i < cloneset[i].nVar(); ++seg_i) {
                        ofs << cloneset[i].getVarGeneStart(seg_i) << "|";
                        ofs << cloneset[i].getVarSeqStart(seg_i) << "|";
                        ofs << cloneset[i].getVarLen(seg_i);
                        CELL_FILL(seg_i, cloneset[i].nVar(), ofs, ',', '\t')
                    }

                    if (gene_segments.is_vdj()) {
                        //
                        // TODO: make so every D alignment will show up in the output file
                        //
                        for (auto seg_i = 0; seg_i < cloneset[i].nDiv(); ++seg_i) {
                            ofs << cloneset[i].getDivGeneStart(seg_i, 0) << "|";
                            ofs << cloneset[i].getDivSeqStart(seg_i, 0) << "|";
                            ofs << cloneset[i].getDivLen(seg_i, 0);
                            CELL_FILL(seg_i, cloneset[i].nDiv(), ofs, ';', '\t')
                        }
                    } else {
                        ofs << "\t";
                    }

                    for (auto seg_i = 0; seg_i < cloneset[i].nJoi(); ++seg_i) {
                        ofs << cloneset[i].getJoiGeneStart(seg_i) << "|";
                        ofs << cloneset[i].getJoiSeqStart(seg_i) << "|";
                        ofs << cloneset[i].getJoiLen(seg_i);
                        CELL_FILL(seg_i, cloneset[i].nJoi(), ofs, ',', std::endl)
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
