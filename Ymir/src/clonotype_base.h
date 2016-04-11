//
// Created by Vadim N. on 11/04/2016.
//

#ifndef YMIR_CLONOTYPE_BASE_H
#define YMIR_CLONOTYPE_BASE_H


#include "alignment.h"
#include "codon_table.h"


namespace ymir {


    class ClonotypeBase {
    public:

        ClonotypeBase(const sequence_t &sequence,
                      SequenceType seq_type,
                      Recombination recomb,
                      bool good)
            : _sequence(sequence),
              _recomb(recomb),
              _seq_type(seq_type),
              _good(good)
        {
        }


        virtual ~ClonotypeBase()
        {
        }


        /**
         * \brief Get the sequence of this Clonotype.
         *
         * Function sequence() returns just stored sequence;
         * nuc_sequence() tries to return the nucleotide sequence and fails in debug mode
         * if the stored sequence is the amino acid one;
         * aa_sequence() always returns amino acid sequence.
         */
        ///@{
        const sequence_t& sequence() const { return _sequence; }

        const sequence_t& nuc_sequence() const {
#ifndef DNDEBUG
            check_and_throw(_seq_type != NUCLEOTIDE, "Clonotype's call to nuc_sequence() is incorrect: wrong sequence type.");
#endif
            return _sequence;
        }

        sequence_t aa_sequence() const {
#ifndef DNDEBUG
            check_and_throw(_seq_type == UNDEF_SEQ_TYPE, "Clonotype's call to aa_sequence() is incorrect: undefined sequence type.");
#endif
            if (_seq_type == NUCLEOTIDE) {
                return translate(_sequence);
            } else {
                return _sequence;
            }
        }
        ///@}


        sequence_t::const_iterator seq_iterator(seq_len_t pos) const { return _sequence.cbegin() + pos; }


        Recombination recombination() const { return _recomb; }


        SequenceType sequence_type() const { return _seq_type; }


        /**
         * \brief Check if clonotype's sequence is coding, noncoding or out-of-frame.
         */
        ///@{
        bool isCoding() const {
            if (_seq_type == NUCLEOTIDE) {
                return !(is_out_of_frame(_sequence) || has_end_codon(_sequence));
            } else {
                return !has_bad_aa_codons(_sequence);
            }
        }

        bool isNoncoding() const {
            if (_seq_type == NUCLEOTIDE) {
                return is_out_of_frame(_sequence) || has_end_codon(_sequence);
            } else {
                return has_bad_aa_codons(_sequence);
            }
        }

        bool isOutOfFrame() const {
            if (_seq_type == NUCLEOTIDE) {
                return is_out_of_frame(_sequence);
            } else {
                return has_oof_aa_codon(_sequence);
            }
        }
        ///@}

        /**
         * \brief Return the string representation of the clonotype.
         */
//        std::string toString() const {
//            std::string res = "";
//            res += this->nuc_sequence() + "\t";
//            res += "V:" + std::to_string(_segments[0]) + "\t";
//            res += "D:" + std::to_string(_segments[2]) + "\t";
//            res += "J:" + std::to_string(_segments[1]);
//            return res;
//        }

        bool is_good() const { return _good; }



    protected:

        sequence_t _sequence; //* CDR3 or full nucleotide or amino acid sequence of a clone. */

        Recombination _recomb;

        SequenceType _seq_type;

        bool _good;


        /**
         *
         */
        ClonotypeBase() {}

    };
}


#endif //YMIR_CLONOTYPE_BASE_H
