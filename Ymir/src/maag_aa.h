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

        friend class MAAGBuilder;
        friend class MAAGForwardBackwardAlgorithm;

    public:

        /**
         * \brief Default constructor.
         */
        MAAGaa()
                : MAAGBase(AMINOACID) //, codons, insertions
        {
            _values.reserve(1);
        }

        /**
         *
         */
        MAAGaa(const MAAGaa &other)
                : MAAGBase(other),
                  _codons(other._codons),
                  _insertions(other._insertions->clone())
        {

        }


        MAAGaa(MAAGaa &&other)
                : MAAGBase(other)
        {
            _codons.swap(other._codons);
            _insertions.swap(other._insertions);
        }


        /**
         *
         */
        virtual ~MAAGaa() { }


        MAAGaa& operator= (const MAAGaa &other) {
            MAAGBase::operator=(other);

            // codons

            // insertions

            return *this;
        }


        MAAGaa& operator=(MAAGaa &&other) {
            MAAGBase::operator=(other);

            // codons

            // insertions

            return *this;
        }


        /**
         * \brief Compute the full generation probability of the amino acid sequence.
         */
        prob_t fullProbability(MAAGComputeProbAction action = SUM_PROBABILITY) const {
#ifndef DNDEBUG
            assert(_insertions);
#endif

        }



    protected:

        CodonMMC _codons;
        unique_ptr<AbstractInsertionModel> _insertions;

    };

}


#endif //YMIR_MAAG_AA_H
