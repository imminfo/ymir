//
// Created by Vadim N. on 24/03/2015.
//

#ifndef _YMIR_CLONOTYPEASSEMBLER_H_
#define _YMIR_CLONOTYPEASSEMBLER_H_


#include "insertionmodel.h"
#include "clonotype.h"


namespace ymir {

    class ClonotypeAssembler {

        struct VJGenerator {
        public:

            VJGenerator() {

            }

            Clonotype generate(ModelParameterVector param_vec) {
//                v_del_num = std::discrete_distribution<eventind_t>(param_vec.get_iterator(param_vec.event_index(VJ_VAR_DEL, vgene, 0)),
//                                                                   param_vec.get_iterator(param_vec.event_index(VJ_VAR_DEL, vgene, 0))
//                                                                   + param_vec.eventFamilySize(VJ_VAR_DEL, vgene - 1))(_rg);
//
//                j_del_num = std::discrete_distribution<eventind_t>(param_vec.get_iterator(param_vec.event_index(VJ_JOI_DEL, jgene, 0)),
//                                                                   param_vec.get_iterator(param_vec.event_index(VJ_JOI_DEL, jgene, 0))
//                                                                   + param_vec.eventFamilySize(VJ_JOI_DEL, jgene - 1))(_rg);


            }

        private:
            std::default_random_engine _rg;
            std::discrete_distribution<segindex_t> _genes;

        };

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
            std::default_random_engine rg;

            if (_param_vec.recombination() == VJ_RECOMB) {
                for (size_t clonotype_i = 0; clonotype_i < count; ++clonotype_i) {
                    vec.push_back(this->generate_vj(rg));
                    if ((clonotype_i + 1) % 50000 == 0) {
                        cout << "Generated " << (clonotype_i + 1) << " sequences." << endl;
                    }
                }
            }
            else {
                const InsertionModel mc_vd(DI_NUCLEOTIDE, _param_vec.get_iterator(_param_vec.event_prob(VDJ_VAR_DIV_INS_NUC, 0, 0)));
                const InsertionModel mc_dj(DI_NUCLEOTIDE, _param_vec.get_iterator(_param_vec.event_prob(VDJ_DIV_JOI_INS_NUC, 0, 0)));

                for (size_t clonotype_i = 0; clonotype_i < count; ++clonotype_i) {

                    //
                    // GENERATION IS HERE
                    //

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


        Clonotype generate_vj(std::default_random_engine &rg) const {
            ClonotypeBuilder builder;
            InsertionModel mc_vj(MONO_NUCLEOTIDE, _param_vec.get_iterator(_param_vec.event_prob(VJ_VAR_JOI_INS_NUC, 0, 0)));

//            std::discrete_distribution<segindex_t> vgenes(_param_vec);
//            vgenes(rg);
//            mc_vj.generate(rg);

            return builder.buildClonotype();
        }

    };

}

#endif //_YMIR_CLONOTYPEASSEMBLER_H_
