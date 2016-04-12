//
// Created by Vadim N. on 09/04/2016.
//

#ifndef YMIR_MAAG_AA_H
#define YMIR_MAAG_AA_H


#include "maag_base.h"


namespace ymir {

    class MAAGaa;


    /**
     * \class MAAGaa
     */
    class MAAGaa : public MAAGBase {
    public:

        /**
         * \brief Default constructor.
         */
        MAAGaa()
                : MAAGBase(AMINOACID) //, codons
        {
            _values.reserve(1);
        }

        /**
         *
         */
        MAAGaa(const MAAGaa &other)
                : MAAGBase(other) //, codons
        {
        }


        MAAGaa(MAAGaa &&other) {
            // codons
        }


        /**
         *
         */
        virtual ~MAAGaa() { }


        MAAGaa& operator= (const MAAGaa &other) {
            MAAGBase::operator=(other);

            // codons

            return *this;
        }


        MAAGaa& operator=(MAAGaa &&other) {
            MAAGBase::operator=(other);

            // codons

            return *this;
        }


        /**
         * \brief Compute the full generation probability of the amino acid sequence.
         */
        prob_t fullProbability(MAAGComputeProbAction action = SUM_PROBABILITY) const {

        }



    protected:


    };

}


#endif //YMIR_MAAG_AA_H
