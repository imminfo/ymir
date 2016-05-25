//
// Created by Vadim N. on 11/04/2016.
//

#ifndef YMIR_CLONOTYPE_NUC_H
#define YMIR_CLONOTYPE_NUC_H


#include "clonotype_base.h"
#include "vdj_alignment_nuc.h"


namespace ymir {


    class ClonotypeNuc : public VDJAlignmentNuc, public ClonotypeBase {

    public:

        static const SequenceType sequence_type = NUCLEOTIDE;


        typedef VDJAlignmentNuc vdj_alignment_t;


        ClonotypeNuc(const sequence_t &sequence,
                     Recombination recomb,
                     const segments_storage_t &segments,
                     const NoGapAlignmentVector &alignments,
                     const n_D_alignments_storage_t &n_D_alignments)
            : ClonotypeBase(sequence,
                            recomb,
                            segments[0] && segments[1] && ((segments[2] && recomb == VDJ_RECOMB) || recomb == VJ_RECOMB)),
              VDJAlignmentNuc(segments,
                              alignments,
                              n_D_alignments)
        {
        }


        ClonotypeNuc(const sequence_t &sequence,
                     Recombination recomb,
                     const VDJAlignmentNuc &alignment)
            : ClonotypeBase(sequence,
                            recomb,
                            alignment.nVar() && alignment.nJoi() && ((alignment.nDiv() && recomb == VDJ_RECOMB) || recomb == VJ_RECOMB)),
              VDJAlignmentNuc(alignment)
        {
        }


        virtual ~ClonotypeNuc()
        {
        }


        sequence_t aa_sequence() const { return translate(_sequence); }


        /**
         * \brief Check if clonotype's sequence is coding, noncoding or out-of-frame.
         */
        ///@{
        bool isCoding() const { return !(is_out_of_frame(_sequence) || has_end_codon(_sequence)); }

        bool isNoncoding() const { return is_out_of_frame(_sequence) || has_end_codon(_sequence); }

        bool isOutOfFrame() const { return is_out_of_frame(_sequence); }
        ///@}


        std::string toString(const VDJRecombinationGenes &genes) const {
            std::string res = _sequence;

            res += "\n";

            for (int i = 0; i < this->nVar(); ++i) {
                for (int j = 0; j < this->getVarSeqStart(i) - 1; ++j) {
                    res += " ";
                }
                for (int j = this->getVarGeneStart(i); j <= this->getVarGeneEnd(i); ++j) {
                    res += genes.V()[this->getVar(i)].sequence[j - 1];
                }
                res += "\n";

                for (int j = 0; j < this->getVarSeqStart(i) - 1; ++j) {
                    res += " ";
                }
                for (int j = this->getVarSeqStart(i); j <= this->getVarSeqEnd(i); ++j) {
                    res += this->isVarMismatch(i, j) ? "1" : "0";
                }
                res += "\n";

                res += " | "
                       + std::to_string(this->getVar(i))
                       + "("
                       + std::to_string(this->getVarSeqStart(i))
                       + ";"
                       + std::to_string(this->getVarSeqEnd(i))
                       + ") ("
                       + std::to_string(this->getVarGeneStart(i))
                       + ";"
                       + std::to_string(this->getVarGeneEnd(i))
                       + ")";
                res += "\n";
            }

//            for (int i = 0; i < this->nJoi(); ++i) {
//                res += "\n";
//            }

            return res;
        }


    protected:

        ClonotypeNuc()
        {
        }

    };


    typedef unique_ptr<ClonotypeNuc> ClonotypeNucPtr;

}


#endif //YMIR_CLONOTYPE_NUC_H
