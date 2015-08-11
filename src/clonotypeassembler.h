//
// Created by Vadim N. on 24/03/2015.
//

#ifndef _YMIR_CLONOTYPEASSEMBLER_H_
#define _YMIR_CLONOTYPEASSEMBLER_H_


#include <chrono>

#include "insertionmodel.h"
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

            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine rg(seed);

            if (_param_vec.recombination() == VJ_RECOMB) {
                std::cout << "Generating sequences:" << std::endl;

                for (size_t clonotype_i = 0; clonotype_i < count; ++clonotype_i) {
                    vec.push_back(this->generate_vj(rg));

                    if ((clonotype_i + 1) % 50000 == 0) {
                        cout << "Generated " << (clonotype_i + 1) << " sequences." << endl;
                    }
                }
            }
            else if (_param_vec.recombination() == VDJ_RECOMB) {
                std::cout << "Generating sequences:" << std::endl;

                const InsertionModel mc_vd(DI_NUCLEOTIDE, _param_vec.get_iterator(_param_vec.event_prob(VDJ_VAR_DIV_INS_NUC, 0, 0)));
                const InsertionModel mc_dj(DI_NUCLEOTIDE, _param_vec.get_iterator(_param_vec.event_prob(VDJ_DIV_JOI_INS_NUC, 0, 0)));

                for (size_t clonotype_i = 0; clonotype_i < count; ++clonotype_i) {
//                    vec.push_back(this->generate_vdj(rg));

                    if ((clonotype_i + 1) % 50000 == 0) {
                        cout << "Generated " << (clonotype_i + 1) << " sequences." << endl;
                    }
                }
            } else {
                std::cout << "Unrecognised recombination type of the input model." << std::endl;
            }

            return Cloneset(vec);
        }

    protected:

        ModelParameterVector _param_vec;
        VDJRecombinationGenes _genes;


        Clonotype generate_vj(std::default_random_engine &rg) const {
            ClonotypeBuilder builder;
            InsertionModel mc_vj(MONO_NUCLEOTIDE, _param_vec.get_iterator(_param_vec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)));

            std::discrete_distribution<event_ind_t> vj_genes(_param_vec.get_iterator(_param_vec.event_index(VJ_VAR_JOI_GEN, 0, 0)),
                                                             _param_vec.get_iterator(_param_vec.event_index(VJ_VAR_JOI_GEN, 0, 0)
                                                                                     + _param_vec.eventClassSize(VJ_VAR_JOI_GEN)));
            event_ind_t v_j = vj_genes(rg);

            seg_index_t vgene = v_j / _param_vec.n_columns(VJ_VAR_JOI_GEN) + 1;

            seg_index_t jgene = v_j - (v_j / _param_vec.n_columns(VJ_VAR_JOI_GEN)) * (_param_vec.n_columns(VJ_VAR_JOI_GEN)) + 1;

            seq_len_t v_del_num = std::discrete_distribution<event_ind_t>(_param_vec.get_iterator(_param_vec.event_index(VJ_VAR_DEL, vgene - 1, 0)),
                                                                         _param_vec.get_iterator(_param_vec.event_index(VJ_VAR_DEL, vgene - 1, 0)
                                                                         + _param_vec.eventFamilySize(VJ_VAR_DEL, vgene - 1)))(rg);
            builder.addValignment(vgene, _genes.V()[vgene].sequence.size() - v_del_num);

            seq_len_t ins_len = std::discrete_distribution<seq_len_t>(_param_vec.get_iterator(_param_vec.event_index(VJ_VAR_JOI_INS_LEN, 0, 0)),
                                                                      _param_vec.get_iterator(_param_vec.event_index(VJ_VAR_JOI_INS_LEN, 0, 0)
                                                                                              + _param_vec.eventClassSize(VJ_VAR_JOI_INS_LEN)))(rg);

            seq_len_t j_del_num = std::discrete_distribution<event_ind_t>(_param_vec.get_iterator(_param_vec.event_index(VJ_JOI_DEL, jgene - 1, 0)),
                                                                         _param_vec.get_iterator(_param_vec.event_index(VJ_JOI_DEL, jgene - 1, 0)
                                                                         + _param_vec.eventFamilySize(VJ_JOI_DEL, jgene - 1)))(rg);
            builder.addJalignment(jgene, _genes.V()[vgene].sequence.size() - v_del_num + ins_len + 1);


//            std::string vgene_del = _genes.V()[vgene].sequence.substr(0, _genes.V()[vgene].sequence.size() - v_del_num);
//            std::string jgene_del = _genes.J()[jgene].sequence.substr(j_del_num);

//            std::string ins_seq = mc_vj.generate(ins_len, rg);
//
//            std::cout << "V+Ji:\t" << (int) v_j << std::endl;
//            std::cout << "V:\t" << _genes.V()[vgene].sequence << std::endl;
//            std::cout << "Vi:\t" << (event_ind_t) vgene << std::endl;
//            std::cout << "#V del:\t" << (int) v_del_num << std::endl;
//            std::cout << "V+del:\t" << vgene_del << std::endl;
//            std::cout << "ins len:\t" << ins_len << std::endl;
//            std::cout << "ins nuc:\t" << ins_seq << std::endl;
//            std::cout << "J:\t" << _genes.J()[jgene].sequence << std::endl;
//            std::cout << "Ji:\t" << (event_ind_t) jgene << std::endl;
//            std::cout << "#J del:\t" << (int) j_del_num << std::endl;
//            std::cout << "J+del:\t" << jgene_del << std::endl;
//
//            std:string seq = vgene_del + "|" + ins_seq + "|" + jgene_del;
//            std::cout << "seq:\t" << seq << std::endl;
//            builder.setSequence(seq);
            builder.setSequence(_genes.V()[vgene].sequence.substr(0, _genes.V()[vgene].sequence.size() - v_del_num)
                                  + mc_vj.generate(ins_len, rg)
                                  + _genes.J()[jgene].sequence.substr(j_del_num));
            builder.setNucleotideSeq();
            return builder.buildClonotype();
        }

    };

}

#endif //_YMIR_CLONOTYPEASSEMBLER_H_
