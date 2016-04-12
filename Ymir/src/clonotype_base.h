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
                      Recombination recomb,
                      bool good)
            : _sequence(sequence),
              _recomb(recomb),
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
         */
        const sequence_t& sequence() const { return _sequence; }


        sequence_t::const_iterator seq_iterator(seq_len_t pos) const { return _sequence.cbegin() + pos; }


        Recombination recombination() const { return _recomb; }

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

        bool _good;


        /**
         *
         */
        ClonotypeBase() {}

    };
}


#endif //YMIR_CLONOTYPE_BASE_H
