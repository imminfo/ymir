//
// Created by Vadim N. on 24/03/2015.
//

#ifndef _YMIR_CLONOTYPEASSEMBLER_H_
#define _YMIR_CLONOTYPEASSEMBLER_H_


#include <random>

#include "clonotype.h"


namespace ymir {

    class ClonotypeAssembler {
    public:

        /**
         *
         */
        ClonotypeAssembler(const ModelParameterVector &param_vec, const VDJRecombinationGenes &genes)
                : _param_vec(param_vec), _genes(genes)
        { }


        /**
         *
         */
        virtual ~ClonotypeAssembler() { }


        /**
         *
         */
        Cloneset generate(size_t count = 1) const {
            std::vector<Clonotype> vec;
            vec.reserve(count);
            ClonotypeBuilder builder;
            string sequence;
            std::default_random_engine rg;

            if (_param_vec.recombination() == VJ_RECOMB) {
                const MarkovChain mc_vj(_param_vec.get_iterator(_param_vec.event_prob(VJ_VAR_JOI_INS_NUC, 0, 0)));
                std::uniform_int_distribution<int> vgenes(1,6);
                vgenes(rg);

                for (size_t clonotype_i = 0; clonotype_i < count; ++clonotype_i) {
                    sequence = "";

                    //
                    // GENERATION IS HERE
                    //

                    vec.push_back(builder.buildClonotype());
                    if ((clonotype_i + 1) % 50000 == 0) {
                        cout << "Generated " << (clonotype_i + 1) << " sequences." << endl;
                    }
                }
            }
            else {
                const MarkovChain mc_vd(_param_vec.get_iterator(_param_vec.event_prob(VDJ_VAR_DIV_INS_NUC, 0, 0)));
                const MarkovChain mc_dj(_param_vec.get_iterator(_param_vec.event_prob(VDJ_DIV_JOI_INS_NUC, 0, 0)));

                for (size_t clonotype_i = 0; clonotype_i < count; ++clonotype_i) {
                    sequence = "";

                    //
                    // GENERATION IS HERE
                    //

                    vec.push_back(builder.buildClonotype());
                    if ((clonotype_i + 1) % 50000 == 0) {
                        cout << "Generated " << (clonotype_i + 1) << " sequences." << endl;
                    }
                }
            }

            return Cloneset(vec);
        }

    protected:

        ModelParameterVector _param_vec;
        VDJRecombinationGenes _genes;

    };

}

#endif //_YMIR_CLONOTYPEASSEMBLER_H_
