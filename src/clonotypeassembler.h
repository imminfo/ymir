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
                        cout << "Generated " << (size_t) (clonotype_i + 1) << "/" << (size_t) count << " sequences." << endl;
                    }
                }
                cout << "Generated " << (size_t) count << "/" << (size_t) count << " sequences." << endl;
            }
            else if (_param_vec.recombination() == VDJ_RECOMB) {
                std::cout << "Generating sequences:" << std::endl;

                for (size_t clonotype_i = 0; clonotype_i < count; ++clonotype_i) {
                    vec.push_back(this->generate_vdj(rg));

                    if ((clonotype_i + 1) % 50000 == 0) {
                        cout << "Generated " << (size_t) (clonotype_i + 1) << "/" << (size_t) count << " sequences." << endl;
                    }
                }
                cout << "Generated " << (size_t) count << "/" << (size_t) count << " sequences." << endl;
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
            builder.addJalignment(jgene, _genes.J()[jgene].sequence.size() - j_del_num + ins_len + 1);

            builder.setSequence(_genes.V()[vgene].sequence.substr(0, _genes.V()[vgene].sequence.size() - v_del_num)
                                  + mc_vj.generate(ins_len, rg)
                                  + _genes.J()[jgene].sequence.substr(j_del_num));
            builder.setNucleotideSeq();
            return builder.buildClonotype();
        }


        Clonotype generate_vdj(std::default_random_engine &rg) const {
            ClonotypeBuilder builder;
            const InsertionModel mc_vd(DI_NUCLEOTIDE, _param_vec.get_iterator(_param_vec.event_index(VDJ_VAR_DIV_INS_NUC, 0, 0)));
            const InsertionModel mc_dj(DI_NUCLEOTIDE, _param_vec.get_iterator(_param_vec.event_index(VDJ_DIV_JOI_INS_NUC, 0, 0)));

            seq_len_t ins_len_vd = std::discrete_distribution<seq_len_t>(_param_vec.get_iterator(_param_vec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 0)),
                                                                         _param_vec.get_iterator(_param_vec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 0)
                                                                                                 + _param_vec.eventClassSize(VDJ_VAR_DIV_INS_LEN)))(rg);
            seq_len_t ins_len_dj = std::discrete_distribution<seq_len_t>(_param_vec.get_iterator(_param_vec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 0)),
                                                                         _param_vec.get_iterator(_param_vec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 0)
                                                                                                 + _param_vec.eventClassSize(VDJ_DIV_JOI_INS_LEN)))(rg);

            seg_index_t vgene = std::discrete_distribution<event_ind_t>(_param_vec.get_iterator(_param_vec.event_index(VDJ_VAR_GEN, 0, 0)),
                                                                        _param_vec.get_iterator(_param_vec.event_index(VDJ_VAR_GEN, 0, 0)
                                                                                                + _param_vec.eventClassSize(VDJ_VAR_GEN)))(rg) + 1;

            std::discrete_distribution<event_ind_t> jd_genes(_param_vec.get_iterator(_param_vec.event_index(VDJ_JOI_DIV_GEN, 0, 0)),
                                                             _param_vec.get_iterator(_param_vec.event_index(VDJ_JOI_DIV_GEN, 0, 0)
                                                                                     + _param_vec.eventClassSize(VDJ_JOI_DIV_GEN)));
            event_ind_t j_d = jd_genes(rg);

            seg_index_t jgene = j_d / _param_vec.n_columns(VDJ_JOI_DIV_GEN) + 1;

            seg_index_t dgene = j_d - (j_d / _param_vec.n_columns(VDJ_JOI_DIV_GEN)) * _param_vec.n_columns(VDJ_JOI_DIV_GEN) + 1;

            seq_len_t v_del_num = std::discrete_distribution<event_ind_t>(_param_vec.get_iterator(_param_vec.event_index(VDJ_VAR_DEL, vgene - 1, 0)),
                                                                          _param_vec.get_iterator(_param_vec.event_index(VDJ_VAR_DEL, vgene - 1, 0)
                                                                                                  + _param_vec.eventFamilySize(VDJ_VAR_DEL, vgene - 1)))(rg);
            std::string vgene_del = _genes.V()[vgene].sequence.substr(0, _genes.V()[vgene].sequence.size() - v_del_num);
            builder.addValignment(vgene, vgene_del.size());

            event_ind_t d_del_index = std::discrete_distribution<event_ind_t>(_param_vec.get_iterator(_param_vec.event_index(VDJ_DIV_DEL, dgene - 1, 0)),
                                                                              _param_vec.get_iterator(_param_vec.event_index(VDJ_DIV_DEL, dgene - 1, 0)
                                                                                                  + _param_vec.eventFamilySize(VDJ_DIV_DEL, dgene - 1)))(rg);
            seq_len_t d5_del_num = d_del_index / _param_vec.n_columns(VDJ_DIV_DEL, 0);
            seq_len_t d3_del_num = d_del_index - (d_del_index / _param_vec.n_columns(VDJ_DIV_DEL, 0)) * _param_vec.n_columns(VDJ_DIV_DEL, 0);
            std::string dgene_del = _genes.D()[dgene].sequence.substr(d5_del_num, _genes.D()[dgene].sequence.size() - d3_del_num - d5_del_num - 1);
            builder.addDalignment(dgene, vgene_del.size() + ins_len_vd + 1,
                                  vgene_del.size() + ins_len_vd + 1 + dgene_del.size(),
                                  d3_del_num + 1,
                                  _genes.D()[dgene].sequence.size() - d3_del_num);

            seq_len_t j_del_num = std::discrete_distribution<event_ind_t>(_param_vec.get_iterator(_param_vec.event_index(VDJ_JOI_DEL, jgene - 1, 0)),
                                                                          _param_vec.get_iterator(_param_vec.event_index(VDJ_JOI_DEL, jgene - 1, 0)
                                                                                                  + _param_vec.eventFamilySize(VDJ_JOI_DEL, jgene - 1)))(rg);
            std::string jgene_del = _genes.J()[jgene].sequence.substr(j_del_num);
            builder.addJalignment(jgene, _genes.J()[jgene].sequence.size() - j_del_num + ins_len_vd + ins_len_dj + dgene_del.size() + 1);
            
            builder.setSequence(vgene_del
                                + mc_vd.generate(ins_len_vd, rg, vgene_del[vgene_del.size() - 1])
                                + dgene_del
                                + mc_dj.generate(ins_len_dj, rg, jgene_del[0], true)
                                + jgene_del);

            builder.setNucleotideSeq();
            return builder.buildClonotype();
        }

    };

}

#endif //_YMIR_CLONOTYPEASSEMBLER_H_
