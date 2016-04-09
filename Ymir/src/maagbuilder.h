/*
 * Ymir <imminfo.github.io/ymir>
 *
 * This file is part of Ymir, a fast C++ tool for computation of assembling
 * probabilities, statistical inference of assembling statistical model
 * and generation of artificial sequences of T-cell receptors data.
 *
 *
 * Copyright 2015 Vadim Nazarov <vdn at mailbox dot com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef _MAAGBUILDER_H_
#define _MAAGBUILDER_H_

#define DEFAULT_SEQ_POSES_RESERVE 300

#define VARIABLE_GENES_MATRIX_INDEX 0
#define VARIABLE_DELETIONS_MATRIX_INDEX 1
#define JOINING_DELETIONS_VJ_MATRIX_INDEX 3
#define JOINING_GENES_VDJ_MATRIX_INDEX 6
#define JOINING_DELETIONS_VDJ_MATRIX_INDEX 5
#define DIVERSITY_GENES_MATRIX_INDEX 3
#define VarJoi_INSERTIONS_MATRIX_INDEX 2
#define VarDiv_INSERTIONS_MATRIX_INDEX 2
#define DivJoi_INSERTIONS_MATRIX_INDEX 4


#include "clonotype.h"
#include "genesegment.h"
#include "insertionmodel.h"
#include "maag.h"
#include "modelparametervector.h"
#include "repertoire.h"


namespace ymir {

    class MAAGBuilder;


    /**
    * \class MAAGBuilder
    */
    class MAAGBuilder : protected MAAG {

    public:

//        static const size_t verbose_step = 10000;


        /**
         * \brief Constructor for the builder from given vector with event probabilities and gene segments.
         */
        MAAGBuilder(const ModelParameterVector &param_vec, const VDJRecombinationGenes &genes)
                : MAAG(),
                  _param_vec(new ModelParameterVector(param_vec)),
                  _genes(new VDJRecombinationGenes(genes)) {
        }


        virtual ~MAAGBuilder() {
        }


        void updateModelParameterVector(const ModelParameterVector &param_vec) {
            *(_param_vec.get()) = param_vec;
        }


        /**
         * \brief Build MAAGs from the given clonotypes.
         *
         * \param clonotype Clonotype from which build the MAAG.
         * \param cloneset Set of clonotypes from which build the repertoires of MAAGs.
         * \param metadata_mode If true than make MAAG with stored event indices.
         * \param aminoacid If true than build MAAGs from the aminoacid sequences.
         *
         * \return Newly constructed MAAG.
         */
        ///@{
        MAAG build(const Clonotype &clonotype,
                   MetadataMode metadata_mode,
                   ErrorMode error_mode,
                   SequenceType seq_type = NUCLEOTIDE) const {
            if (clonotype.is_good()) {
                ProbMMC probs;
                EventIndMMC events;
                ErrMMC errors;
                vector<seq_len_t> seq_poses;
                seq_poses.reserve(DEFAULT_SEQ_POSES_RESERVE);

                auto resize_size = 0, e_resize_size = 0;
                switch (clonotype.recombination()) {
                    case VJ_RECOMB:
                        resize_size = VJ_CHAIN_SIZE;
                        e_resize_size = 2;
                        break;

                    case VDJ_RECOMB:
                        resize_size = VDJ_CHAIN_SIZE;
                        e_resize_size = 3;
                        break;

                    default:
#ifndef DNDEBUG
                        check_and_throw(false, "MAAGBuilder: unknown recombination type.");
#endif
                }

                probs.resize(resize_size);
                if (metadata_mode) {
                    events.resize(resize_size);
                }
                if (error_mode) {
                    errors.resize(e_resize_size);
                }

                this->buildVariable(clonotype, probs, events, errors, seq_poses, metadata_mode, error_mode);
                this->buildJoining(clonotype, probs, events, errors, seq_poses, metadata_mode, error_mode);
                if (clonotype.recombination() == VJ_RECOMB) {
                    this->buildVJinsertions(clonotype, probs, events, seq_poses, metadata_mode, error_mode);
                } else if (clonotype.recombination() == VDJ_RECOMB) {
                    this->buildDiversity(clonotype, probs, events, errors, seq_poses, metadata_mode, error_mode);
                    this->buildVDinsertions(clonotype, probs, events, seq_poses, metadata_mode, error_mode);
                    this->buildDJinsertions(clonotype, probs, events, seq_poses, metadata_mode, error_mode);
                }

//              VERY OLD VERSION
//            if (error_mode && metadata_mode) {
//
//                // TODO: deal with D deletions and insertions null matrices
//                // - if the D deletions matrix contains only zeros, then remove this matrix
//                // - if insertions matrices have columns / rows with only zero (!) events (!), then remove it
//                // and remove the corresponding deletions rows / columns from neighbour matrices.
//                //
//                if (clonotype.recombination() == VDJ_RECOMB) {
//
//                }
//
//                unique_ptr<seq_len_t[]> seq_poses_arr(new seq_len_t[seq_poses.size()]);
//                copy(seq_poses.begin(), seq_poses.end(), seq_poses_arr.get());
//                return MAAG(probs, events, errors, clonotype.sequence(), seq_poses_arr, seq_poses.size(), seq_type);
//            } else if (metadata_mode) {
//                unique_ptr<seq_len_t[]> seq_poses_arr(new seq_len_t[seq_poses.size()]);
//                copy(seq_poses.begin(), seq_poses.end(), seq_poses_arr.get());
//                return MAAG(probs, events, clonotype.sequence(), seq_poses_arr, seq_poses.size(), seq_type);
//            } else if (error_mode) {
//                return MAAG(probs, errors);
//            } else {
//                return MAAG(probs);
//            }

//                std::cout << clonotype.toString() << std::endl;


                // OLD VERSION
                probs.finish();

                MAAG maag;
                maag._recomb = clonotype.recombination();
                maag._seq_type = clonotype.sequence_type();
                maag.swap(probs);
                if (error_mode) {
                    errors.finish();

                    maag._errors.reset(new ErrMMC());
                    maag._errors->swap(errors);
                }
                if (metadata_mode) {
                    events.finish();

                    unique_ptr<seq_len_t[]> seq_poses_arr(new seq_len_t[seq_poses.size()]);
                    copy(seq_poses.begin(), seq_poses.end(), seq_poses_arr.get());
                    maag._sequence.reset(new sequence_t(clonotype.sequence()));
                    maag._seq_poses.swap(seq_poses_arr);
                    maag._n_poses = seq_poses.size();

                    maag._events.reset(new EventIndMMC());
                    maag._events->swap(events);
                }
                return maag;


                // NEW VERSION
//                probs->finish();
//
//                MAAG maag;
//                maag._recomb = clonotype.recombination();
//                maag._seq_type = clonotype.sequence_type();
//                maag.swap(*probs);
//                if (error_mode) {
//                    errors->finish();
//
//                    maag._errors.swap(errors);
//                }
//                if (metadata_mode) {
//                    events->finish();
//
//                    unique_ptr<seq_len_t[]> seq_poses_arr(new seq_len_t[seq_poses.size()]);
//                    copy(seq_poses.begin(), seq_poses.end(), seq_poses_arr.get());
//                    maag._sequence.reset(new sequence_t(clonotype.sequence()));
//                    maag._seq_poses.swap(seq_poses_arr);
//                    maag._n_poses = seq_poses.size();
//
//                    maag._events.swap(events);
//                }
//                return maag;

        } else {
            return MAAG();
        }

    }

        MAAGRepertoire build(const ClonesetView &cloneset,
                             MetadataMode metadata_mode,
                             ErrorMode error_mode,
                             SequenceType seq_type = NUCLEOTIDE,
                             bool verbose = true) const
        {
            size_t verbose_step;

            if (verbose) {
                std::cout << "Building " << (size_t) cloneset.size() << " MAAGs..." << std::endl;
                verbose_step = cloneset.size() / 10;
            }

            MAAGRepertoire res;
            res.resize(cloneset.size());
#ifdef USE_OMP
#if OMP_THREADS == -1
            #pragma omp parallel for 
#else
            #pragma omp parallel for  num_threads(OMP_THREADS)
#endif
#endif
            for (size_t i = 0; i < cloneset.size(); ++i) {
                res[i] = this->build(cloneset[i], metadata_mode, error_mode, seq_type);

#ifndef USE_OMP
                if (verbose && (i+1) % verbose_step == 0 && (i+1) != cloneset.size()) {
                    cout << "[" << (int) ((100*(i+1)) / cloneset.size()) << "%] "<< "Built " << (int) (i+1) << " / " << (int) (cloneset.size()) << " MAAGs." << endl;
                }
#endif
            }

            if (verbose) {
                cout << "[100%] Built " << (int) (cloneset.size()) << " MAAGs." << endl;
            }

            return res;
        }
        ///@}


        /**
         * \brief Compute generation probabilities without building the full information about MAAGs.
         *
         * \param clonotype
         * \param cloneset
         * \param aminoacid If true then compute amino acid generation probabilities.
         * \param action What action to perform for computing generation probabilties - either sum all probabilities
         * or choose the max one.
         */
        ///@{
        prob_t buildAndCompute(const Clonotype &clonotype,
                               ErrorMode error_mode,
                               SequenceType seq_type = NUCLEOTIDE,
                               MAAGComputeProbAction action = SUM_PROBABILITY) const {
            return this->build(clonotype, NO_METADATA, error_mode, seq_type).fullProbability(action);
        }

        vector<prob_t> buildAndCompute(const ClonesetView &cloneset,
                                       ErrorMode error_mode,
                                       SequenceType seq_type = NUCLEOTIDE,
                                       MAAGComputeProbAction action = SUM_PROBABILITY,
                                       bool verbose = true) const
        {
            vector<prob_t> res;
            res.reserve(cloneset.size());
            size_t verbose_step;

            if (verbose) {
                std::cout << "Computing assembling probabilities on " << (size_t) cloneset.size() << " clonotypes." << std::endl;
                verbose_step = cloneset.size() / 10;
            }

#ifdef USE_OMP
#if OMP_THREADS == -1
            #pragma omp parallel for 
#else
            #pragma omp parallel for  num_threads(OMP_THREADS)
#endif
#endif
            for (size_t i = 0; i < cloneset.size(); ++i) {
                res.push_back(buildAndCompute(cloneset[i], error_mode, seq_type, action));

#ifndef USE_OMP
                if (verbose && (i+1) % verbose_step == 0 && (i+1) != cloneset.size()) {
                    std::cout << "[" << (int) ((100*(i+1)) / cloneset.size()) << "%] " << "Computed " << (int) (i+1) << " / " << (size_t) cloneset.size() << " assembling probabilities." << std::endl;
                }
#endif
            }

            if (verbose) {
                std::cout << "[100%] Computed " << (size_t) cloneset.size() << " assembling probabilities." << std::endl;
            }

            return res;
        }
        ///@}


        /**
         * \brief Replace event probabilities in the given MAAGs if they have stored event indices.
         *
         * \param maag MAAG with an event indices matrix.
         * \param repertoire Repertoire with MAAGs with event indices matrices.
         */
        ///@{
        void updateEventProbabilities(MAAG *maag) const {
            if (maag->has_events()) {
                vector<seq_len_t> seq_poses_vec(maag->_seq_poses.get(), maag->_seq_poses.get() + maag->_n_poses);
                for (int node_i = 0; node_i < maag->chainSize(); ++node_i) {
                    // either rebuild all insertions
                    if (maag->is_vj() && node_i == VarJoi_INSERTIONS_MATRIX_INDEX) {
                        MonoNucInsertionModel im(_param_vec->get_iterator(_param_vec->event_index(VJ_VAR_JOI_INS_NUC, 0, 0))); // TODO: add errors here?

                        seq_len_t v_vertices = maag->nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                                j_vertices = maag->nodeRows(JOINING_DELETIONS_VJ_MATRIX_INDEX);

                        this->buildInsertions(maag->sequence(),
                                              *maag,
                                              *maag->_events,
                                              seq_poses_vec,
                                              VarJoi_INSERTIONS_MATRIX_INDEX,
                                              _param_vec->event_index(VJ_VAR_JOI_INS_LEN, 0, 0),
                                              _param_vec->max_VJ_ins_len(),
                                              NO_METADATA,
                                              maag->has_errors(),
                                              0,
                                              v_vertices - 1,
                                              v_vertices,
                                              v_vertices + j_vertices - 1,
                                              im,
                                              false);

                    } else if (maag->is_vdj() && node_i == VarDiv_INSERTIONS_MATRIX_INDEX) {
                        DiNucInsertionModel im(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC, 0, 0)));  // TODO: add errors here?

                        seq_len_t v_vertices = maag->nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                                d3_vertices = maag->nodeRows(DIVERSITY_GENES_MATRIX_INDEX);

                        this->buildInsertions(maag->sequence(),
                                              *maag,
                                              *maag->_events,
                                              seq_poses_vec,
                                              VarDiv_INSERTIONS_MATRIX_INDEX,
                                              _param_vec->event_index(VDJ_VAR_DIV_INS_LEN, 0, 0),
                                              _param_vec->max_VD_ins_len(),
                                              NO_METADATA,
                                              maag->has_errors(),
                                              0,
                                              v_vertices - 1,
                                              v_vertices,
                                              v_vertices + d3_vertices - 1,
                                              im,
                                              false);

                    } else if (maag->is_vdj() && node_i == DivJoi_INSERTIONS_MATRIX_INDEX) {
                        DiNucInsertionModel im(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC, 0, 0)));

                        seq_len_t v_vertices = maag->nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                                d3_vertices = maag->nodeRows(DIVERSITY_GENES_MATRIX_INDEX),
                                d5_vertices = maag->nodeColumns(DIVERSITY_GENES_MATRIX_INDEX),
                                j_vertices = maag->nodeRows(JOINING_DELETIONS_VDJ_MATRIX_INDEX);

                        this->buildInsertions(maag->sequence(),
                                              *maag,
                                              *maag->_events,
                                              seq_poses_vec,
                                              DivJoi_INSERTIONS_MATRIX_INDEX,
                                              _param_vec->event_index(VDJ_DIV_JOI_INS_LEN, 0, 0),
                                              _param_vec->max_DJ_ins_len(),
                                              NO_METADATA,
                                              maag->has_errors(),
                                              v_vertices + d3_vertices,
                                              v_vertices + d3_vertices + d5_vertices - 1,
                                              v_vertices + d3_vertices + d5_vertices,
                                              v_vertices + d3_vertices + d5_vertices + j_vertices - 1,
                                              im,
                                              true);

                    } else {
                        // or just replace all event probabilities with the new ones
                        if (!maag->has_errors()) {
                            for (int mat_i = 0; mat_i < maag->nodeSize(node_i); ++mat_i) {
                                for (int row_i = 0; row_i < maag->nodeRows(node_i); ++row_i) {
                                    for (int col_i = 0; col_i < maag->nodeColumns(node_i); ++col_i) {
                                        (*maag)(node_i, mat_i, row_i, col_i) =
                                                (*_param_vec)[maag->event_index(node_i, mat_i, row_i, col_i)];
                                    }
                                }
                            }
                        } else {
                            int err_node_i = 0;
                            if (maag->is_vj()) {
                                err_node_i = node_i == JOINING_DELETIONS_VJ_MATRIX_INDEX ? 1 : 0;
                            } else {
                                if (node_i != VARIABLE_DELETIONS_MATRIX_INDEX) {
                                    if (node_i == JOINING_DELETIONS_VDJ_MATRIX_INDEX) {
                                        err_node_i = 2;
                                    } else {
                                        err_node_i = 1;
                                    }
                                }
                            }

                            if (node_i != 0 && (maag->is_vdj() && node_i != JOINING_GENES_VDJ_MATRIX_INDEX)) {
                                for (int mat_i = 0; mat_i < maag->nodeSize(node_i); ++mat_i) {
                                    for (int row_i = 0; row_i < maag->nodeRows(node_i); ++row_i) {
                                        for (int col_i = 0; col_i < maag->nodeColumns(node_i); ++col_i) {
                                            if (maag->errors(err_node_i, mat_i, row_i, col_i)) {
                                                (*maag)(node_i, mat_i, row_i, col_i) =
                                                        (*_param_vec)[maag->event_index(node_i, mat_i, row_i, col_i)]
                                                        * _param_vec->error_prob() * maag->errors(err_node_i, mat_i, row_i, col_i);
                                            } else {
                                                (*maag)(node_i, mat_i, row_i, col_i) = (*_param_vec)[maag->event_index(node_i, mat_i, row_i, col_i)];
                                            }
                                        }
                                    }
                                }
                            } else {
                                for (int mat_i = 0; mat_i < maag->nodeSize(node_i); ++mat_i) {
                                    for (int row_i = 0; row_i < maag->nodeRows(node_i); ++row_i) {
                                        for (int col_i = 0; col_i < maag->nodeColumns(node_i); ++col_i) {
                                            (*maag)(node_i, mat_i, row_i, col_i) = (*_param_vec)[maag->event_index(node_i, mat_i, row_i, col_i)];
                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }

        void updateEventProbabilities(MAAGRepertoire *repertoire, bool verbose = true) const {
            size_t verbose_step;

            if (verbose) {
                std::cout << "Updating " << (size_t) repertoire->size() << " MAAGs..." << std::endl;
                verbose_step = repertoire->size() / 10;
            }

#ifdef USE_OMP
#if OMP_THREADS == -1
            #pragma omp parallel for 
#else
            #pragma omp parallel for  num_threads(OMP_THREADS)
#endif
#endif
            for (size_t i = 0; i < repertoire->size(); ++i) {
                this->updateEventProbabilities(&(*repertoire)[i]);

#ifndef USE_OMP
                if (verbose && (i+1) % verbose_step == 0 && (i+1) != repertoire->size()) {
                    cout << "[" << (int) ((100*(i+1)) / repertoire->size()) << "%] " << "Updated " << (size_t) (i+1) << " / " << (size_t) repertoire->size() << " MAAGs." << endl;
                }
#endif
            }

            if (verbose) {
                cout << "[100%] Updated " << (int) (repertoire->size()) << " MAAGs." << endl;
            }
        }
        ///@}


    protected:

        unique_ptr<ModelParameterVector> _param_vec;  // or just copy it?
        unique_ptr<VDJRecombinationGenes> _genes; // copy this too?


        /**
         * \brief Private default constructor.
         */
        MAAGBuilder() : _param_vec(nullptr), _genes(nullptr) {}


        /**
        * \brief Build probability and events matrices for Variable gene segments.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param metadata_mode Boolean if build should be full.
        */
        void buildVariable(const Clonotype &clonotype,
                           ProbMMC &probs,
                           EventIndMMC &events,
                           ErrMMC &errors,
                           vector<seq_len_t> &seq_poses,
                           MetadataMode metadata_mode,
                           ErrorMode error_mode) const
        {
            // find max V alignment
            seq_len_t len = 0;
            seg_index_t v_num = clonotype.nVar(), j_num = clonotype.nJoi();
            for (int v_index = 0; v_index < v_num; ++v_index) {
                len = std::max(len, clonotype.getVarLen(v_index));
            }

            // compute V deletions
            seq_len_t v_len, v_gene, v_start, v_end;

            probs.initNode(VARIABLE_DELETIONS_MATRIX_INDEX, v_num, 1, len + 1);
            if (metadata_mode) {
                events.initNode(VARIABLE_DELETIONS_MATRIX_INDEX, v_num, 1, len + 1);
            }
            if (error_mode) {
                errors.initNode(0, v_num, 1, len + 1);
            }

            if (clonotype.recombination() == VJ_RECOMB) {
                probs.initNode(VARIABLE_GENES_MATRIX_INDEX, 1, v_num, j_num);
                if (metadata_mode) {
                    events.initNode(VARIABLE_GENES_MATRIX_INDEX, 1, v_num, j_num);
                }
            } else if (clonotype.recombination() == VDJ_RECOMB) {
                probs.initNode(VARIABLE_GENES_MATRIX_INDEX, v_num, 1, 1);
                if (metadata_mode) {
                    events.initNode(VARIABLE_GENES_MATRIX_INDEX, v_num, 1, 1);
                }
            }

            EventClass V_DEL = clonotype.recombination() == VJ_RECOMB ? VJ_VAR_DEL : VDJ_VAR_DEL;
            for (seg_index_t v_index = 0; v_index < v_num; ++v_index) {

                // probability of choosing this V gene segment
                v_gene = clonotype.getVar(v_index);

                if (clonotype.recombination() == VJ_RECOMB) {
                    for (seg_index_t j_index = 0; j_index < j_num; ++j_index) {
                        probs(VARIABLE_GENES_MATRIX_INDEX, 0, v_index, j_index)
                                = _param_vec->event_prob(VJ_VAR_JOI_GEN, 0, v_gene - 1, clonotype.getJoi(j_index) - 1);
                    }
                } else if (clonotype.recombination() == VDJ_RECOMB) {
                    probs(VARIABLE_GENES_MATRIX_INDEX, v_index, 0, 0) = _param_vec->event_prob(VDJ_VAR_GEN, 0, v_gene - 1); // probability of choosing this V gene segment
                }

                // V deletions
                v_len = _genes->V()[v_gene].sequence.size();
                v_start = clonotype.getVarGeneStart(v_index);
                v_end = clonotype.getVarGeneEnd(v_index);

                for (seq_len_t i = 0; i < v_end - v_start + 2; ++i) {
                    probs(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = _param_vec->event_prob(V_DEL, v_gene - 1, (1 + v_len) - (v_start + i)); // probability of deletions
                }

                if (metadata_mode) {
                    if (clonotype.recombination() == VJ_RECOMB) {
                        // probability of choosing this V gene segment
                        for (seg_index_t j_index = 0; j_index < j_num; ++j_index) {
                            events(VARIABLE_GENES_MATRIX_INDEX, 0, v_index, j_index)
                                    = _param_vec->event_index(VJ_VAR_JOI_GEN, 0, v_gene - 1, clonotype.getJoi(j_index) - 1);
                        }
                    } else if (clonotype.recombination() == VDJ_RECOMB) {
                        events(VARIABLE_GENES_MATRIX_INDEX, v_index, 0, 0) = _param_vec->event_index(VDJ_VAR_GEN, 0, v_gene - 1);
                    }

                    for (seq_len_t i = 0; i < v_end - v_start + 2; ++i) {
                        events(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) = _param_vec->event_index(V_DEL, v_gene - 1, (1 + v_len) - (v_start + i));
                    }
                }

                if (error_mode) {
                    errors(0, v_index, 0, 0) = 0;
                    for (seq_len_t i = 1; i <= clonotype.getVarLen(v_index); ++i) {
                        errors(0, v_index, 0, i) = errors(0, v_index, 0, i-1) + clonotype.isVarMismatch(v_index, i);
                        if (errors(0, v_index, 0, i)) {
                            probs(VARIABLE_DELETIONS_MATRIX_INDEX, v_index, 0, i) *= errors(0, v_index, 0, i) * _param_vec->error_prob();
                        }
                    }
                }
            }

            for (seq_len_t i = 0; i <= len; ++i) { seq_poses.push_back(i); }
        }


        /**
        * \brief Build probability and events matrices for Joining gene segments.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param metadata_mode Boolean if build should be full.
        */
        void buildJoining(const Clonotype &clonotype,
                          ProbMMC &probs,
                          EventIndMMC &events,
                          ErrMMC &errors,
                          vector<seq_len_t> &seq_poses,
                          MetadataMode metadata_mode,
                          ErrorMode error_mode) const
        {
            int J_index_dels = JOINING_DELETIONS_VJ_MATRIX_INDEX,
                    J_index_genes = JOINING_GENES_VDJ_MATRIX_INDEX;
            if (clonotype.recombination() == VDJ_RECOMB) {
                J_index_dels = JOINING_DELETIONS_VDJ_MATRIX_INDEX;
            }

            // find max J alignment
            seg_index_t j_num = clonotype.nJoi();
            seq_len_t len = 0, seq_global_start_pos = (seq_len_t) -1;
            for (int j_index = 0; j_index < j_num; ++j_index) {
                len = std::max(len, clonotype.getJoiLen(j_index));
                seq_global_start_pos = std::min(seq_global_start_pos, clonotype.getJoiSeqStart(j_index));
            }

            // add J deletions nodes
            probs.initNode(J_index_dels, j_num, len + 1, 1);
            if (metadata_mode) {
                events.initNode(J_index_dels, j_num, len + 1, 1);
            }
            if (error_mode) {
                errors.initNode(errors.chainSize() - 1, j_num, len + 1, 1);
            }

            // add J or J-D gene nodes
            if (clonotype.recombination() == VDJ_RECOMB) {
                probs.initNode(J_index_genes, 1, j_num, clonotype.nDiv());
                if (metadata_mode) {
                    events.initNode(J_index_genes, 1, j_num, clonotype.nDiv());
                }
            }

            // compute J deletions
            seq_len_t j_len, j_gene, j_start, j_end, shift;

            EventClass J_DEL = clonotype.recombination() == VJ_RECOMB ? VJ_JOI_DEL : VDJ_JOI_DEL;
            for (seg_index_t j_index = 0; j_index < j_num; ++j_index) {
                // probability of choosing the J segment
                j_gene = clonotype.getJoi(j_index);
                j_len = _genes->J()[j_gene].sequence.size();

                if (clonotype.recombination() == VDJ_RECOMB) {
                    for (seg_index_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
                        probs(J_index_genes, 0, j_index, d_index)
                                = _param_vec->event_prob(VDJ_JOI_DIV_GEN, 0, j_gene - 1, clonotype.getDiv(d_index) - 1); // probability of choosing this J gene segment with other D genes
                    }
                }

                // J deletions
                j_start = clonotype.getJoiGeneStart(j_index);
                j_end = clonotype.getJoiGeneEnd(j_index);
                shift = clonotype.getJoiSeqStart(j_index) - seq_global_start_pos;

                for (seq_len_t i = 0; i < clonotype.getJoiLen(j_index) + 1; ++i) {
                     probs(J_index_dels, j_index, i + shift, 0) = _param_vec->event_prob(J_DEL, j_gene - 1, j_start + i - 1); // probability of deletions
                }

                if (metadata_mode) {
                    if (clonotype.recombination() == VDJ_RECOMB) {
                        for (seg_index_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
                            events(J_index_genes, 0, j_index, d_index)
                                    = _param_vec->event_index(VDJ_JOI_DIV_GEN, 0, j_gene - 1, clonotype.getDiv(d_index) - 1); // probability of choosing this J gene segment with other D genes
                        }
                    }

                    for (seq_len_t i = 0; i < clonotype.getJoiLen(j_index) + 1; ++i) {
                        events(J_index_dels, j_index, i + shift, 0) = _param_vec->event_index(J_DEL, j_gene - 1, j_start + i - 1);
                    }
                }

                if (error_mode) {
                    errors(errors.chainSize() - 1, j_index, len, 0) = 0;
                    for (seq_len_t i = 1; i <= clonotype.getJoiLen(j_index); ++i) {
                        errors(errors.chainSize() - 1, j_index, len - i, 0)
                                = errors(errors.chainSize() - 1, j_index, len + 1 - i, 0)
                                  + clonotype.isJoiMismatch(j_index, clonotype.getJoiLen(j_index) + 1 - i);
                        if (errors(errors.chainSize() - 1, j_index, len - i, 0)) {
                            probs(J_index_dels, j_index, len - i, 0) *= errors(errors.chainSize() - 1, j_index, len - i, 0) * _param_vec->error_prob();
                        }
                    }
                }
            }

            for (seq_len_t i = clonotype.sequence().size() - len + 1; i <= clonotype.sequence().size() + 1; ++i) {
                seq_poses.push_back(i);
            }
        }


        /**
        * \brief Build probability and events matrices for Diversity gene segments.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param metadata_mode Boolean if build should be full.
        */
        void buildDiversity(const Clonotype &clonotype,
                            ProbMMC &probs,
                            EventIndMMC &events,
                            ErrMMC &errors,
                            vector<seq_len_t> &seq_poses,
                            MetadataMode metadata_mode,
                            ErrorMode error_mode) const
        {
            seq_len_t min_D_len;

            // vector seq_start -> 0 means no such index in the matrix, 1 otherwise.
            seq_len_t seq_arr_size = clonotype.sequence().size() + 1;
            unique_ptr<seq_len_t[]> seq_row(new seq_len_t[seq_arr_size]);
            std::fill(seq_row.get(), seq_row.get() + seq_arr_size, 0);

            // vector seq_end -> 0 means no such index in the matrix, 1 otherwise.
            unique_ptr<seq_len_t[]> seq_col(new seq_len_t[seq_arr_size]);
            std::fill(seq_col.get(), seq_col.get() + seq_arr_size, 0);

            for (seg_index_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
                min_D_len = _param_vec->D_min_len(clonotype.getDiv(d_index));

                for (seg_index_t j = 0; j < clonotype.numDivAlignments(d_index); ++j) {
                    seq_len_t d_seq_start = clonotype.getDivSeqStart(d_index, j),
                              d_seq_end = clonotype.getDivSeqEnd(d_index, j);

                    // yes-yes, I know that it could be done more efficiently. But I don't want to.
                    for (seq_len_t i = d_seq_start; i <= d_seq_end - min_D_len + 1; ++i) {
                        seq_row[i] = 1;
                    }

                    for (seq_len_t i = d_seq_start + min_D_len - (seq_len_t) 1; i <= d_seq_end; ++i) {
                        seq_col[i] = 1;
                    }
                }
            }

            // make new vector seq_start -> 1based index in rows of Ddel matrices
            seq_len_t seq_ind = 1;
            for (seq_len_t i = 0; i < seq_arr_size; ++i) {
                if (seq_row[i]) {
                    seq_row[i] = seq_ind;
                    ++seq_ind;
                }
            }

            // make new vector seq_end -> 1based index in columns of Ddel matrices
            seq_ind = 1;
            for (seq_len_t i = 0; i < seq_arr_size; ++i) {
                if (seq_col[i]) {
                    seq_col[i] = seq_ind;
                    ++seq_ind;
                }
            }

            // find indices of D alignments and use only them to reduce memory usage.
            seq_len_t last_max_seq_start = 0, last_max_seq_end = 0, seq_row_ind = 0, seq_col_ind = 0;

            seq_len_t seq_row_nonzeros = 0, seq_col_nonzeros = 0;
            for (seq_len_t i = 0; i < seq_arr_size; ++i) {
                seq_row_nonzeros += seq_row[i] != 0;
                seq_col_nonzeros += seq_col[i] != 0;
            }

            probs.initNode(DIVERSITY_GENES_MATRIX_INDEX, clonotype.nDiv(), seq_row_nonzeros, seq_col_nonzeros);
            if (metadata_mode) {
                events.initNode(DIVERSITY_GENES_MATRIX_INDEX, clonotype.nDiv(), seq_row_nonzeros, seq_col_nonzeros);
            }
            if (error_mode) {
                errors.initNode(1, clonotype.nDiv(), seq_row_nonzeros, seq_col_nonzeros);
            }


            seg_index_t d_index, d_gene;
            seq_len_t d_len;
            seq_len_t d_seq_start, d_seq_end, d_gene_start, d_gene_end;

            for (seg_index_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
                d_gene = clonotype.getDiv(d_index);
                d_len = _genes->D()[d_gene].sequence.size();
                min_D_len = _param_vec->D_min_len(d_gene);

                // for each aligned Div segment get all possible smaller alignments and add them to the matrix.
                for (seg_index_t j = 0; j < clonotype.numDivAlignments(d_index); ++j) {
                    d_seq_start = clonotype.getDivSeqStart(d_index, j);
                    d_seq_end = clonotype.getDivSeqEnd(d_index, j);
                    d_gene_start = clonotype.getDivGeneStart(d_index, j);
                    d_gene_end = clonotype.getDivGeneEnd(d_index, j);

                    for (seq_len_t left_pos = d_seq_start; left_pos <= d_seq_end - min_D_len + 1; ++left_pos) {
                        for (seq_len_t right_pos = left_pos + min_D_len - 1; right_pos <= d_seq_end; ++right_pos) {
                            probs(DIVERSITY_GENES_MATRIX_INDEX, d_index, seq_row[left_pos] - 1, seq_col[right_pos] - 1)
                                    = _param_vec->event_prob(VDJ_DIV_DEL,
                                                             d_gene - 1,
                                                             d_gene_start + left_pos - d_seq_start,
                                                             d_len - (d_gene_end - (d_seq_end - right_pos)));
                            if (metadata_mode) {
                                events(DIVERSITY_GENES_MATRIX_INDEX, d_index, seq_row[left_pos] - 1, seq_col[right_pos] - 1)
                                        = _param_vec->event_index(VDJ_DIV_DEL,
                                                                  d_gene - 1,
                                                                  d_gene_start + left_pos - d_seq_start,
                                                                  d_len - (d_gene_end - (d_seq_end - right_pos)));
                            }

                            if (error_mode) {
                                errors(1, d_index, seq_row[left_pos] - 1, seq_col[right_pos] - 1) =
                                        clonotype.numDivMismatches(d_index, j,
                                                                   d_gene_start + left_pos - d_seq_start,
                                                                   d_gene_end - (d_seq_end - right_pos));

                                if (errors(1, d_index, seq_row[left_pos] - 1, seq_col[right_pos] - 1)) {
                                    probs(DIVERSITY_GENES_MATRIX_INDEX, d_index, seq_row[left_pos] - 1, seq_col[right_pos] - 1)
                                            *= _param_vec->error_prob() * errors(1, d_index, seq_row[left_pos] - 1, seq_col[right_pos] - 1);
                                }
                            }
                        }
                    }
                }
            }

            // insert D3 and D5 positions
            vector<seq_len_t> D35_poses;
            D35_poses.reserve(seq_row_nonzeros + seq_col_nonzeros + 2);
            for (seq_len_t i = 1; i < seq_arr_size; ++i) { if (seq_row[i]) { D35_poses.push_back(i); } }
            for (seq_len_t i = 1; i < seq_arr_size; ++i) { if (seq_col[i]) { D35_poses.push_back(i); } }

            // Note! insert diversity gene seq poses BEFORE joining gene seq poses
            seq_poses.reserve(seq_poses.size() + D35_poses.size() + 2);  // +2 -> just in case (:
            seq_poses.insert(seq_poses.begin() + probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX), D35_poses.begin(), D35_poses.end());
        }


        /**
        * \brief Build probability and events matrices for Variable-Joining gene segments insertions.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param metadata_mode Boolean if build should be full.
        */
        void buildVJinsertions(const Clonotype &clonotype,
                               ProbMMC &probs,
                               EventIndMMC &events,
                               const vector<seq_len_t> &seq_poses,
                               MetadataMode metadata_mode,
                               ErrorMode error_mode) const
        {
            MonoNucInsertionModel mc(_param_vec->get_iterator(_param_vec->event_index(VJ_VAR_JOI_INS_NUC, 0, 0)), error_mode ? _param_vec->error_prob() : 0);

            seq_len_t v_vertices = probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                    j_vertices = probs.nodeRows(JOINING_DELETIONS_VJ_MATRIX_INDEX);

            probs.initNode(VarJoi_INSERTIONS_MATRIX_INDEX, 1, v_vertices, j_vertices);

            if (metadata_mode) {
                events.initNode(VarJoi_INSERTIONS_MATRIX_INDEX, 1, v_vertices, j_vertices);
            }

            this->buildInsertions(clonotype.sequence(),
                                  probs,
                                  events,
                                  seq_poses,
                                  VarJoi_INSERTIONS_MATRIX_INDEX,
                                  _param_vec->event_index(VJ_VAR_JOI_INS_LEN, 0, 0),
                                  _param_vec->max_VJ_ins_len(),
                                  metadata_mode,
                                  error_mode,
                                  0,
                                  v_vertices - 1,
                                  v_vertices,
                                  v_vertices + j_vertices - 1,
                                  mc,
                                  false);
        }


        /**
        * \brief Build probability and events matrices for Variable-Diversity gene segments insertions.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param metadata_mode Boolean if build should be full.
        */
        void buildVDinsertions(const Clonotype &clonotype,
                               ProbMMC &probs,
                               EventIndMMC &events,
                               const vector<seq_len_t> &seq_poses,
                               MetadataMode metadata_mode,
                               ErrorMode error_mode) const
        {
            DiNucInsertionModel mc(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC, 0, 0)), error_mode ? _param_vec->error_prob() : 0);

            seq_len_t v_vertices = probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                    d3_vertices = probs.nodeRows(DIVERSITY_GENES_MATRIX_INDEX);

            probs.initNode(VarDiv_INSERTIONS_MATRIX_INDEX, 1, v_vertices, d3_vertices);

            if (metadata_mode) {
                events.initNode(VarDiv_INSERTIONS_MATRIX_INDEX, 1, v_vertices, d3_vertices);
            }

            this->buildInsertions(clonotype.sequence(),
                                  probs,
                                  events,
                                  seq_poses,
                                  VarDiv_INSERTIONS_MATRIX_INDEX,
                                  _param_vec->event_index(VDJ_VAR_DIV_INS_LEN, 0, 0),
                                  _param_vec->max_VD_ins_len(),
                                  metadata_mode,
                                  error_mode,
                                  0,
                                  v_vertices - 1,
                                  v_vertices,
                                  v_vertices + d3_vertices - 1,
                                  mc,
                                  false);
        }


        /**
        * \brief Build probability and events matrices for Diversity-Joining gene segments insertions.
        *
        * \param clonotype Clonotype that used for building the graph.
        * \param probs Multi-Matrix Chain with event probabilities.
        * \param events Multi-Matrix Chain with event indices.
        * \param seq_poses Vector of positions.
        * \param metadata_mode Boolean if build should be full.
        */
        void buildDJinsertions(const Clonotype &clonotype,
                               ProbMMC &probs,
                               EventIndMMC &events,
                               const vector<seq_len_t> &seq_poses,
                               MetadataMode metadata_mode,
                               ErrorMode error_mode) const
        {
            DiNucInsertionModel mc(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC, 0, 0)), error_mode ? _param_vec->error_prob() : 0);

            seq_len_t v_vertices = probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                    d3_vertices = probs.nodeRows(DIVERSITY_GENES_MATRIX_INDEX),
                    d5_vertices = probs.nodeColumns(DIVERSITY_GENES_MATRIX_INDEX),
                    j_vertices = probs.nodeRows(JOINING_DELETIONS_VDJ_MATRIX_INDEX);

            probs.initNode(DivJoi_INSERTIONS_MATRIX_INDEX, 1, d5_vertices, j_vertices);

            if (metadata_mode) {
                events.initNode(DivJoi_INSERTIONS_MATRIX_INDEX, 1, d5_vertices, j_vertices);
            }

            this->buildInsertions(clonotype.sequence(),
                                  probs,
                                  events,
                                  seq_poses,
                                  DivJoi_INSERTIONS_MATRIX_INDEX,
                                  _param_vec->event_index(VDJ_DIV_JOI_INS_LEN, 0, 0),
                                  _param_vec->max_DJ_ins_len(),
                                  metadata_mode,
                                  error_mode,
                                  v_vertices + d3_vertices,
                                  v_vertices + d3_vertices + d5_vertices - 1,
                                  v_vertices + d3_vertices + d5_vertices,
                                  v_vertices + d3_vertices + d5_vertices + j_vertices - 1,
                                  mc,
                                  true);
        }


        /**
         * \brief General function for building insertions.
         *
         * \param clonotype Clonotype that used for building the graph.
         * \param probs Multi-Matrix Chain with event probabilities.
         * \param events Multi-Matrix Chain with event indices.
         * \param seq_poses Vector of positions.
         * \param ins_node_index Node of the event / prob matrix with insertions.
         * \param null_insertion Event index of insertions of length zero.
         * \param max_size Maximum size of the length of insertions.
         * \param metadata_mode Boolean if build should be full.
         * \param left_vertices_start Starting index in seq_poses for the vertices in the left matrix.
         * \param left_vertices_end Ending index in seq_poses for the vertices in the left matrix.
         * \param right_vertices_start Starting index in seq_poses for the vertices in the right matrix.
         * \param right_vertices_end Ending index in seq_poses for the vertices in the right matrix.
         * \param mc Insertion model that uses for generation of N nucleotides.
        */
        void buildInsertions(const string &sequence,
                             ProbMMC &probs,
                             EventIndMMC &events,
                             const vector<seq_len_t> &seq_poses,
                             ProbMMC::node_ind_t ins_node_index,
                             event_ind_t null_insertion,
                             seq_len_t max_size,
                             MetadataMode metadata_mode,
                             bool error_mode,
                             seq_len_t left_vertices_start,
                             seq_len_t left_vertices_end,
                             seq_len_t right_vertices_start,
                             seq_len_t right_vertices_end,
                             const AbstractInsertionModel& mc,
                             bool reversed = false) const
        {
            int insertion_len;
            bool good_insertion;
            char last_char = NULL_CHAR;

            for (size_t left_vertex_i = left_vertices_start; left_vertex_i <= left_vertices_end; ++left_vertex_i) {
                for (size_t right_vertex_i = right_vertices_start; right_vertex_i <= right_vertices_end; ++right_vertex_i) {
                    insertion_len = seq_poses[right_vertex_i] - seq_poses[left_vertex_i] - 1;
                    good_insertion = (insertion_len >= 0) && (insertion_len <= max_size);
                    if (good_insertion) {
                        last_char = NULL_CHAR;
                        if (!reversed) {
                            if (seq_poses[left_vertex_i] != 0) {
                                last_char = sequence[seq_poses[left_vertex_i] - 1];
                            }

                            probs(ins_node_index, 0, left_vertex_i - left_vertices_start, right_vertex_i - right_vertices_start)
                                    = mc.nucProbability(sequence.cbegin() + seq_poses[left_vertex_i],
                                                        insertion_len,
                                                        last_char,
                                                        error_mode)
                                      * (*_param_vec)[null_insertion + insertion_len];
                        } else {
                            if (seq_poses[right_vertex_i] != sequence.size() + 1) {
                                last_char = sequence[seq_poses[right_vertex_i] - 1];
                            }

                            probs(ins_node_index, 0, left_vertex_i - left_vertices_start, right_vertex_i - right_vertices_start)
                                    = mc.nucProbability(sequence.crbegin() + (sequence.size() - seq_poses[right_vertex_i] + 1),
                                                        insertion_len,
                                                        last_char,
                                                        error_mode)
                                      * (*_param_vec)[null_insertion + insertion_len];
                        }

                        if (metadata_mode) {
                            events(ins_node_index, 0, left_vertex_i - left_vertices_start, right_vertex_i - right_vertices_start)
                                    = null_insertion + insertion_len;
                        }
                    }
                }
            }
        }


        // build<VJ_RECOMB, SAVE_METADATA, NO_ERR>
        // put specific values (lambda) to a specific MMC
        // function for finding max V alignment
        // function for finding max J alignment
        // function for shrinking D alignment matrices, find non-zero positions, etc.
        // make seq_row, seq_col, seq_start, seq_ind vector (one function?)
        // build[Mono|Di]NucInsertions <InsertionModel, SequenceType, MetadataMode>
        // general functions for assigning values (event probs / event inds) to MMC of some type.

        /*
            buildVarGenesAndDels
            buildDivDels
            buildJoiDels
            buildJoiDivGenes
        */

    };
}

#endif //_MAAGBUILDER_H_
