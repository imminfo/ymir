//
// Created by Vadim N. on 12/04/2016.
//

#ifndef YMIR_MAAG_BUILDER_H
#define YMIR_MAAG_BUILDER_H


#include "genesegment.h"
#include "insertion_model.h"
#include "maag.h"
#include "modelparametervector.h"
#include "repertoire.h"


namespace ymir {

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


    class MAAGBuilder {
    public:

        /**
         * \brief Constructor for the builder from given vector with event probabilities and gene segments.
         */
        MAAGBuilder(const ModelParameterVector &param_vec, const VDJRecombinationGenes &genes)
            : _param_vec(new ModelParameterVector(param_vec)),
              _genes(new VDJRecombinationGenes(genes))
        {
             this->update_insertion_models();
        }


        void updateModelParameterVector(const ModelParameterVector &param_vec) {
            *(_param_vec.get()) = param_vec;
             this->update_insertion_models();
        }


        //
        // Nucleotide sequences
        //

        /**
         * \brief Build MAAGs from the given clonotypes.
         *
         * \param clonotype Clonotype from which build the MAAG.
         * \param cloneset Set of clonotypes from which build the repertoires of MAAGs.
         * \param metadata_mode If true than make MAAG with stored event indices.
         *
         * \return Newly constructed MAAG.
         */
        ///@{
        MAAGnuc build(const ClonotypeNuc &clonotype, MetadataMode metadata_mode, ErrorMode error_mode) const;

        MAAGNucRepertoire build(const ClonesetViewNuc &cloneset, MetadataMode metadata_mode, ErrorMode error_mode, bool verbose = true) const;
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
        prob_t buildAndCompute(const ClonotypeNuc &clonotype, ErrorMode error_mode, MAAGComputeProbAction action = SUM_PROBABILITY) const;

        vector<prob_t> buildAndCompute(const ClonesetViewNuc &cloneset, ErrorMode error_mode, MAAGComputeProbAction action = SUM_PROBABILITY, bool verbose = true) const;
        ///@}


        /**
         * \brief Replace event probabilities in the given MAAGs if they have stored event indices.
         *
         * \param maag MAAG with an event indices matrix.
         * \param repertoire Repertoire with MAAGs with event indices matrices.
         */
        ///@{
        void updateEventProbabilities(MAAGnuc *maag) const;

        void updateEventProbabilities(MAAGNucRepertoire *repertoire, bool verbose = true) const;
        ///@}


        //
        // Amino acid sequences
        //

        /**
         *
         */
        ///@{
        MAAGaa build(const ClonotypeAA &clonotype) const;

        MAAGAARepertoire build(const ClonesetViewAA &cloneset, bool verbose = true) const;
        ///@}


        /**
         *
         */
        ///@{
        prob_t buildAndCompute(const ClonotypeAA &clonotype) const;

        vector<prob_t> buildAndCompute(const ClonesetViewAA &cloneset, bool verbose = true) const;
        ///@}

    private:

        unique_ptr<ModelParameterVector> _param_vec;
        unique_ptr<VDJRecombinationGenes> _genes;
        unique_ptr<MonoNucInsertionModel> _vj_ins;
        unique_ptr<DiNucInsertionModel> _vd_ins, _dj_ins;


        /**
         * \brief Private default constructor.
         */
        MAAGBuilder() : _param_vec(nullptr), _genes(nullptr) { }


        void update_insertion_models() {
            if (_param_vec) {
                assert(_param_vec->recombination() == VJ_RECOMB || _param_vec->recombination() == VDJ_RECOMB);
                if (_param_vec->recombination() == VJ_RECOMB) {
                    _vj_ins.reset(new MonoNucInsertionModel(_param_vec->get_iterator(_param_vec->event_index(VJ_VAR_JOI_INS_NUC, 0, 0)), _param_vec->error_prob()));
                } else if (_param_vec->recombination() == VDJ_RECOMB) {
                    _vd_ins.reset(new DiNucInsertionModel(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC, 0, 0)), _param_vec->error_prob()));
                    _dj_ins.reset(new DiNucInsertionModel(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC, 0, 0)), _param_vec->error_prob()));
                } else {
                    std::cerr << "Wrong recombination" << std::endl;
                }
            }
        }


        /**
         * \brief Build probability and events matrices for Variable gene segments.
         *
         * \param clonotype Clonotype that used for building the graph.
         * \param probs Multi-Matrix Chain with event probabilities.
         * \param events Multi-Matrix Chain with event indices.
         * \param seq_poses Vector of positions.
         * \param metadata_mode Boolean if build should be full.
         */
        void buildNucVariable(const ClonotypeNuc &clonotype,
                              ProbMMC &probs,
                              EventIndMMC &events,
                              ErrMMC &errors,
                              vector<seq_len_t> &seq_poses,
                              MetadataMode metadata_mode,
                              ErrorMode error_mode) const;


        /**
         * \brief Build probability and events matrices for Joining gene segments.
         *
         * \param clonotype Clonotype that used for building the graph.
         * \param probs Multi-Matrix Chain with event probabilities.
         * \param events Multi-Matrix Chain with event indices.
         * \param seq_poses Vector of positions.
         * \param metadata_mode Boolean if build should be full.
         */
        void buildNucJoining(const ClonotypeNuc &clonotype,
                             ProbMMC &probs,
                             EventIndMMC &events,
                             ErrMMC &errors,
                             vector<seq_len_t> &seq_poses,
                             MetadataMode metadata_mode,
                             ErrorMode error_mode) const;


        /**
         * \brief Build probability and events matrices for Diversity gene segments.
         *
         * \param clonotype Clonotype that used for building the graph.
         * \param probs Multi-Matrix Chain with event probabilities.
         * \param events Multi-Matrix Chain with event indices.
         * \param seq_poses Vector of positions.
         * \param metadata_mode Boolean if build should be full.
         */
        void buildNucDiversity(const ClonotypeNuc &clonotype,
                               ProbMMC &probs,
                               EventIndMMC &events,
                               ErrMMC &errors,
                               vector<seq_len_t> &seq_poses,
                               MetadataMode metadata_mode,
                               ErrorMode error_mode) const;


        /**
         * \brief Build probability and events matrices for Variable-Joining gene segments insertions.
         *
         * \param clonotype Clonotype that used for building the graph.
         * \param probs Multi-Matrix Chain with event probabilities.
         * \param events Multi-Matrix Chain with event indices.
         * \param seq_poses Vector of positions.
         * \param metadata_mode Boolean if build should be full.
         */
        void buildNucVJinsertions(const ClonotypeNuc &clonotype,
                                  ProbMMC &probs,
                                  EventIndMMC &events,
                                  const vector<seq_len_t> &seq_poses,
                                  MetadataMode metadata_mode,
                                  ErrorMode error_mode) const;


        /**
         * \brief Build probability and events matrices for Variable-Diversity gene segments insertions.
         *
         * \param clonotype Clonotype that used for building the graph.
         * \param probs Multi-Matrix Chain with event probabilities.
         * \param events Multi-Matrix Chain with event indices.
         * \param seq_poses Vector of positions.
         * \param metadata_mode Boolean if build should be full.
         */
        void buildNucVDinsertions(const ClonotypeNuc &clonotype,
                                  ProbMMC &probs,
                                  EventIndMMC &events,
                                  const vector<seq_len_t> &seq_poses,
                                  MetadataMode metadata_mode,
                                  ErrorMode error_mode) const;


        /**
         * \brief Build probability and events matrices for Diversity-Joining gene segments insertions.
         *
         * \param clonotype Clonotype that used for building the graph.
         * \param probs Multi-Matrix Chain with event probabilities.
         * \param events Multi-Matrix Chain with event indices.
         * \param seq_poses Vector of positions.
         * \param metadata_mode Boolean if build should be full.
         */
        void buildNucDJinsertions(const ClonotypeNuc &clonotype,
                                  ProbMMC &probs,
                                  EventIndMMC &events,
                                  const vector<seq_len_t> &seq_poses,
                                  MetadataMode metadata_mode,
                                  ErrorMode error_mode) const;


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
        void buildNucInsertions(const string &sequence,
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
                                bool reversed = false) const;


        /**
         *
         */
        void buildAAVariable(const ClonotypeAA &clonotype,
                             ProbMMC &probs,
                             CodonMMC &codons,
                             vector<seq_len_t> &seq_poses) const;


        /**
         *
         */
        void buildAAJoining(const ClonotypeAA &clonotype,
                            ProbMMC &probs,
                            CodonMMC &codons,
                            vector<seq_len_t> &seq_poses) const;


        /**
         *
         */
        void buildAADiversity(const ClonotypeAA &clonotype,
                              ProbMMC &probs,
                              CodonMMC &codons,
                              vector<seq_len_t> &seq_poses) const;

    };


    MAAGnuc MAAGBuilder::build(const ClonotypeNuc &clonotype, MetadataMode metadata_mode, ErrorMode error_mode) const {
//        assert(clonotype.is_good());

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

            this->buildNucVariable(clonotype, probs, events, errors, seq_poses, metadata_mode, error_mode);
            this->buildNucJoining(clonotype, probs, events, errors, seq_poses, metadata_mode, error_mode);
            if (clonotype.recombination() == VJ_RECOMB) {
                this->buildNucVJinsertions(clonotype, probs, events, seq_poses, metadata_mode, error_mode);
            } else if (clonotype.recombination() == VDJ_RECOMB) {
                this->buildNucDiversity(clonotype, probs, events, errors, seq_poses, metadata_mode, error_mode);
                this->buildNucVDinsertions(clonotype, probs, events, seq_poses, metadata_mode, error_mode);
                this->buildNucDJinsertions(clonotype, probs, events, seq_poses, metadata_mode, error_mode);
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

            MAAGnuc maag;
            maag._recomb = clonotype.recombination();
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
            return MAAGnuc();
        }
    }


    MAAGNucRepertoire MAAGBuilder::build(const ClonesetViewNuc &cloneset, MetadataMode metadata_mode,
                                         ErrorMode error_mode, bool verbose) const {
        size_t verbose_step;

        if (verbose) {
            std::cout << "Building " << (size_t) cloneset.size() << " MAAGs..." << std::endl;
            verbose_step = cloneset.size() / 10;
        }

        MAAGNucRepertoire res;
        res.resize(cloneset.size());

        std::chrono::system_clock::time_point tp1 = std::chrono::system_clock::now();
#ifdef USE_OMP
        #pragma omp parallel for
#endif
        for (size_t i = 0; i < cloneset.size(); ++i) {
            res[i] = this->build(cloneset[i], metadata_mode, error_mode);

#ifndef USE_OMP
            if (verbose && (i+1) % verbose_step == 0 && (i+1) != cloneset.size()) {
                cout << "[" << (int) ((100*(i+1)) / cloneset.size()) << "%] "
                     << "Build " << (int) (i+1) << " / " << (int) (cloneset.size()) << " MAAGs. "
                     << print_time(tp1, cloneset.size(), i) << endl;
            }
#endif
        }

        if (verbose) {
            cout << "[100%] Built " << (int) (cloneset.size()) << " MAAGs in " << time_diff_now(tp1) << endl;
        }

        return res;
    }


    prob_t MAAGBuilder::buildAndCompute(const ClonotypeNuc &clonotype, ErrorMode error_mode,
                                        MAAGComputeProbAction action) const {
        return this->build(clonotype, NO_METADATA, error_mode).fullProbability(action);
    }


    vector<prob_t> MAAGBuilder::buildAndCompute(const ClonesetViewNuc &cloneset, ErrorMode error_mode,
                                                MAAGComputeProbAction action, bool verbose) const {
        vector<prob_t> res;
        res.reserve(cloneset.size());
        size_t verbose_step;

        if (verbose) {
            std::cout << "Computing assembling probabilities on " << (size_t) cloneset.size() << " clonotypes (nuc)." << std::endl;
            verbose_step = cloneset.size() / 10;
        }

        std::chrono::system_clock::time_point tp1 = std::chrono::system_clock::now();
#ifdef USE_OMP
        #pragma omp parallel for
#endif
        for (size_t i = 0; i < cloneset.size(); ++i) {
            res.push_back(buildAndCompute(cloneset[i], error_mode, action));

#ifndef USE_OMP
            if (verbose && (i+1) % verbose_step == 0 && (i+1) != cloneset.size()) {
                cout << "[" << (int) ((100*(i+1)) / cloneset.size()) << "%] "
                     << "Computed " << (int) (i+1) << " / " << (int) (cloneset.size()) << " assembling probabilities. "
                     << print_time(tp1, cloneset.size(), i) << endl;
            }
#endif
        }

        if (verbose) {
            std::cout << "[100%] Computed " << (size_t) cloneset.size() << " assembling probabilities in " << time_diff_now(tp1) << std::endl;
        }

        return res;
    }


    void MAAGBuilder::updateEventProbabilities(MAAGnuc *maag) const {
        if (maag->has_events()) {
            vector<seq_len_t> seq_poses_vec(maag->_seq_poses.get(), maag->_seq_poses.get() + maag->_n_poses);
            for (int node_i = 0; node_i < maag->chainSize(); ++node_i) {
                // either rebuild all insertions
                if (maag->is_vj() && node_i == VarJoi_INSERTIONS_MATRIX_INDEX) {
                    seq_len_t v_vertices = maag->nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                            j_vertices = maag->nodeRows(JOINING_DELETIONS_VJ_MATRIX_INDEX);

                    this->buildNucInsertions(maag->sequence(),
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
                                          *_vj_ins,
                                          false);

                } else if (maag->is_vdj() && node_i == VarDiv_INSERTIONS_MATRIX_INDEX) {
                    seq_len_t v_vertices = maag->nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                            d3_vertices = maag->nodeRows(DIVERSITY_GENES_MATRIX_INDEX);

                    this->buildNucInsertions(maag->sequence(),
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
                                          *_vd_ins,
                                          false);

                } else if (maag->is_vdj() && node_i == DivJoi_INSERTIONS_MATRIX_INDEX) {
                    seq_len_t v_vertices = maag->nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                            d3_vertices = maag->nodeRows(DIVERSITY_GENES_MATRIX_INDEX),
                            d5_vertices = maag->nodeColumns(DIVERSITY_GENES_MATRIX_INDEX),
                            j_vertices = maag->nodeRows(JOINING_DELETIONS_VDJ_MATRIX_INDEX);

                    this->buildNucInsertions(maag->sequence(),
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
                                          *_dj_ins,
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


    void MAAGBuilder::updateEventProbabilities(MAAGNucRepertoire *repertoire, bool verbose) const {
        size_t verbose_step;

        if (verbose) {
            std::cout << "Updating " << (size_t) repertoire->size() << " MAAGs..." << std::endl;
            verbose_step = repertoire->size() / 10;
        }

        std::chrono::system_clock::time_point tp1;
        tp1 = std::chrono::system_clock::now();
#ifdef USE_OMP
        #pragma omp parallel for
#endif
        for (size_t i = 0; i < repertoire->size(); ++i) {
            this->updateEventProbabilities(&(*repertoire)[i]);

#ifndef USE_OMP
            if (verbose && (i+1) % verbose_step == 0 && (i+1) != repertoire->size()) {
                cout << "[" << (int) ((100*(i+1)) / repertoire->size()) << "%] "
                     << "Updated " << (int) (i+1) << " / " << (int) (repertoire->size()) << " MAAGs. "
                     << print_time(tp1, repertoire->size(), i) << endl;
            }
#endif
        }

        if (verbose) {
            std::cout << "[100%] Updated " << (int) (repertoire->size()) << " MAAGs in " << time_diff_now(tp1) << std::endl;
        }
    }


    MAAGaa MAAGBuilder::build(const ClonotypeAA &clonotype) const {
        assert(clonotype.is_good());

        if (clonotype.is_good()) {
            ProbMMC probs;
            CodonMMC codons;
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
                    e_resize_size = 4; // 2 matrices for the start and the end of D segments
                    break;

                default:
#ifndef DNDEBUG
                    check_and_throw(true, "MAAGBuilder: unknown recombination type.");
#endif
            }

            probs.resize(resize_size);
            codons.resize(e_resize_size);

            MAAGaa maag;

            this->buildAAVariable(clonotype, probs, codons, seq_poses);
            this->buildAAJoining(clonotype, probs, codons, seq_poses);
            if (clonotype.recombination() == VJ_RECOMB) {
                probs.initNode(VarJoi_INSERTIONS_MATRIX_INDEX,
                               1,
                               probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                               probs.nodeRows(JOINING_DELETIONS_VJ_MATRIX_INDEX));

                maag._insertions.reset(_vj_ins->clone());
            }
            else if (clonotype.recombination() == VDJ_RECOMB) {
                this->buildAADiversity(clonotype, probs, codons, seq_poses);

                probs.initNode(VarDiv_INSERTIONS_MATRIX_INDEX,
                               2,
                               probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                               probs.nodeRows(DIVERSITY_GENES_MATRIX_INDEX));
                probs.initNode(DivJoi_INSERTIONS_MATRIX_INDEX,
                               4,
                               probs.nodeColumns(DIVERSITY_GENES_MATRIX_INDEX),
                               probs.nodeRows(JOINING_DELETIONS_VDJ_MATRIX_INDEX));

                maag._insertions.reset(_vd_ins->clone());
                maag._insertions_rev.reset(_dj_ins->clone());
            }

            probs.finish();

            maag._recomb = clonotype.recombination();
            maag.swap(probs);
            maag._codons.swap(codons);
            unique_ptr<seq_len_t[]> seq_poses_arr(new seq_len_t[seq_poses.size()]);
            std::copy(seq_poses.begin(), seq_poses.end(), seq_poses_arr.get());
            maag._sequence.reset(new sequence_t(clonotype.sequence()));
            maag._seq_poses.swap(seq_poses_arr);
            maag._n_poses = seq_poses.size();

            if (clonotype.recombination() == VJ_RECOMB) {
                maag._ins_start = _param_vec->get_iterator(_param_vec->event_index(VJ_VAR_JOI_INS_LEN, 0, 0));
                maag._max_ins_len = _param_vec->max_VJ_ins_len();
            } else if (clonotype.recombination() == VDJ_RECOMB) {
                maag._ins_start = _param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_LEN, 0, 0));
                maag._ins_start_rev = _param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_LEN, 0, 0));
                maag._max_ins_len = _param_vec->max_VD_ins_len();
                maag._max_ins_len_rev = _param_vec->max_DJ_ins_len();
            }

            return maag;
        } else {
            return MAAGaa();
        }
    }


    MAAGAARepertoire MAAGBuilder::build(const ClonesetViewAA &cloneset, bool verbose) const {
        size_t verbose_step;

        if (verbose) {
            std::cout << "Building " << (size_t) cloneset.size() << " MAAGs..." << std::endl;
            verbose_step = cloneset.size() / 10;
        }

        MAAGAARepertoire res;
        res.resize(cloneset.size());

        std::chrono::system_clock::time_point tp1 = std::chrono::system_clock::now();
#ifdef USE_OMP
        #pragma omp parallel for
#endif
        for (size_t i = 0; i < cloneset.size(); ++i) {
            res[i] = this->build(cloneset[i]);

#ifndef USE_OMP
            if (verbose && (i+1) % verbose_step == 0 && (i+1) != cloneset.size()) {
                cout << "[" << (int) ((100*(i+1)) / cloneset.size()) << "%] "
                     << "Build " << (int) (i+1) << " / " << (int) (cloneset.size()) << " MAAGs. "
                     << print_time(tp1, cloneset.size(), i) << endl;
            }
#endif
        }

        if (verbose) {
            std::cout << "[100%] Built " << (int) (cloneset.size()) << " MAAGs in " << time_diff_now(tp1) << std::endl;
        }

        return res;
    }


    prob_t MAAGBuilder::buildAndCompute(const ClonotypeAA &clonotype) const {
//        std::cout << (int) clonotype.nVar() << std::endl;
//        std::cout << (int) clonotype.nJoi() << std::endl;
//        auto res = this->build(clonotype);
//        return res.fullProbability();
        return this->build(clonotype).fullProbability();
    }


    vector<prob_t> MAAGBuilder::buildAndCompute(const ClonesetViewAA &cloneset, bool verbose) const {
        vector<prob_t> res;
        res.reserve(cloneset.size());
        size_t verbose_step;

        if (verbose) {
            std::cout << "Computing assembling probabilities on " << (size_t) cloneset.size() << " clonotypes (aa)." << std::endl;
            verbose_step = cloneset.size() / 10;
        }

        std::chrono::system_clock::time_point tp1, tp2;
        tp1 = std::chrono::system_clock::now();

#ifdef USE_OMP
        #pragma omp parallel for
#endif
        for (size_t i = 0; i < cloneset.size(); ++i) {
//            std::cout << "i:" << i << std::endl;
//            std::cout << cloneset[i].sequence() << std::endl;
//            std::cout << cloneset[i].recombination() << std::endl;
            res.push_back(buildAndCompute(cloneset[i]));

#ifndef USE_OMP
            if (verbose && (i+1) % verbose_step == 0 && (i+1) != cloneset.size()) {
                tp2 = std::chrono::system_clock::now();
                cout << "[" << (int) ((100*(i+1)) / cloneset.size()) << "%] "
                     << "Computed " << (int) (i+1) << " / " << (int) (cloneset.size()) << " assembling probabilities. "
                     << print_time(tp1, cloneset.size(), i) << endl;
            }
#endif
        }

        if (verbose) {
            std::cout << "[100%] Computed " << (size_t) cloneset.size() << " assembling probabilities in " << time_diff_now(tp1) << std::endl;
        }

        return res;
    }


    void MAAGBuilder::buildNucVariable(const ClonotypeNuc &clonotype, ProbMMC &probs, EventIndMMC &events, ErrMMC &errors,
                                       vector<seq_len_t> &seq_poses, MetadataMode metadata_mode,
                                       ErrorMode error_mode) const {
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


    void MAAGBuilder::buildNucJoining(const ClonotypeNuc &clonotype, ProbMMC &probs, EventIndMMC &events, ErrMMC &errors,
                                      vector<seq_len_t> &seq_poses, MetadataMode metadata_mode,
                                      ErrorMode error_mode) const {
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


    void MAAGBuilder::buildNucDiversity(const ClonotypeNuc &clonotype, ProbMMC &probs, EventIndMMC &events, ErrMMC &errors,
                                        vector<seq_len_t> &seq_poses, MetadataMode metadata_mode,
                                        ErrorMode error_mode) const {
        seq_len_t min_D_len;

        seq_len_t seq_arr_size = clonotype.sequence().size() + 1; // TODO: +2 ??? Because 1-based with the range [0, size+1]
        // vector seq_row -> 0 means no such index in the matrix (row-wise), 1 otherwise
        std::vector<seq_len_t> seq_row(seq_arr_size, 0), seq_row_num = seq_row;

        // vector seq_col -> 0 means no such index in the matrix (column-wise), 1 otherwise.
        std::vector<seq_len_t> seq_col(seq_arr_size, 0);

        // shifts for each pair of positions (i,j) (in case if (i,j) is duplicated)
        std::vector<seq_len_t> shift_mat(seq_arr_size * seq_arr_size);

        // optimization for faster search for the max element in rows
        unique_ptr<bool[]> any_row_element(new bool[seq_arr_size]);

        for (seg_index_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
            min_D_len = _param_vec->D_min_len(clonotype.getDiv(d_index));
            std::fill(shift_mat.begin(), shift_mat.end(), 0);
            std::fill(any_row_element.get(), any_row_element.get() + seq_arr_size, false);

            for (seg_index_t d_align = 0; d_align < clonotype.numDivAlignments(d_index); ++d_align) {
                seq_len_t d_seq_start = clonotype.getDivSeqStart(d_index, d_align),
                        d_seq_end = clonotype.getDivSeqEnd(d_index, d_align);

                for (seq_len_t i = d_seq_start; i <= d_seq_end - min_D_len + 1; ++i) {
                    any_row_element[i] = true;

                    for (seq_len_t j = d_seq_start + min_D_len - static_cast<seq_len_t>(1); j <= d_seq_end; ++j) {
                        ++shift_mat[i*seq_arr_size + j];
                    }
                }

                for (seq_len_t j = d_seq_start + min_D_len - static_cast<seq_len_t>(1); j <= d_seq_end; ++j) {
                    seq_col[j] = 1;
                }
            }

            // for each row get max element, i.e., max number of repeats of a specific row index
            for (seq_len_t i = 0; i < seq_arr_size; ++i) {
                if (any_row_element[i]) {
                    seq_row[i] = std::max(seq_row[i],
                                          *std::max_element(shift_mat.begin() + i*seq_arr_size, shift_mat.begin() + (i+1)*seq_arr_size));
                }
            }
        }

        seq_len_t seq_row_nonzeros = 0, seq_ind = 0;
        std::vector<seq_len_t> seq_row_indices, seq_col_indices;
        seq_row_indices.reserve(seq_arr_size);
        seq_col_indices.reserve(seq_arr_size);

        // make new vector seq_start -> 1based index in rows of Ddel matrices
        for (seq_len_t i = 0; i < seq_arr_size; ++i) {
            if (seq_row[i]) {
                seq_row_indices.push_back(i);

                seq_row_num[i] = seq_row[i];
                seq_row_nonzeros += seq_row[i];

                seq_row[i] = seq_ind;
                seq_ind += seq_row_num[i];
            }
        }

        // make new vector seq_end -> 1based index in columns of Ddel matrices
        seq_ind = 0;
        for (seq_len_t i = 0; i < seq_arr_size; ++i) {
            if (seq_col[i]) {
                seq_col_indices.push_back(i);

                seq_col[i] = seq_ind;
                ++seq_ind;
            }
        }
        seq_len_t seq_col_nonzeros = seq_ind;

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

            std::fill(shift_mat.begin(), shift_mat.end(), 0);

            // for each aligned Div segment get all possible smaller alignments and add them to the matrix.
            for (seg_index_t j = 0; j < clonotype.numDivAlignments(d_index); ++j) {

                d_seq_start = clonotype.getDivSeqStart(d_index, j);
                d_seq_end = clonotype.getDivSeqEnd(d_index, j);
                d_gene_start = clonotype.getDivGeneStart(d_index, j);
                d_gene_end = clonotype.getDivGeneEnd(d_index, j);

                for (seq_len_t left_pos = d_seq_start; left_pos <= d_seq_end - min_D_len + 1; ++left_pos) {
                    for (seq_len_t right_pos = left_pos + min_D_len - 1; right_pos <= d_seq_end; ++right_pos) {
                        probs(DIVERSITY_GENES_MATRIX_INDEX,
                              d_index,
                              seq_row[left_pos] + shift_mat[seq_arr_size*left_pos + right_pos],
                              seq_col[right_pos])
                                =
                                _param_vec->event_prob(VDJ_DIV_DEL,
                                                       d_gene - 1,
                                                       d_gene_start + left_pos - d_seq_start,
                                                       d_len - (d_gene_end - (d_seq_end - right_pos)));
                        if (metadata_mode) {
                            events(DIVERSITY_GENES_MATRIX_INDEX,
                                   d_index,
                                   seq_row[left_pos] + shift_mat[seq_arr_size*left_pos + right_pos],
                                   seq_col[right_pos])
                                    =
                                    _param_vec->event_index(VDJ_DIV_DEL,
                                                            d_gene - 1,
                                                            d_gene_start + left_pos - d_seq_start,
                                                            d_len - (d_gene_end - (d_seq_end - right_pos)));
                        }

                        if (error_mode) {
                            errors(1,
                                   d_index,
                                   seq_row[left_pos] + shift_mat[seq_arr_size*left_pos + right_pos],
                                   seq_col[right_pos])
                                    =
                                    clonotype.numDivMismatches(d_index,
                                                               j,
                                                               d_gene_start + left_pos - d_seq_start,
                                                               d_gene_end - (d_seq_end - right_pos));

                            if (errors(1, d_index, seq_row[left_pos], seq_col[right_pos])) {
                                probs(DIVERSITY_GENES_MATRIX_INDEX,
                                      d_index,
                                      seq_row[left_pos] + shift_mat[seq_arr_size*left_pos + right_pos],
                                      seq_col[right_pos])
                                        *=
                                        _param_vec->error_prob()
                                        * errors(1, d_index, seq_row[left_pos], seq_col[right_pos]);
                            }
                        }

                        ++shift_mat[seq_arr_size*left_pos + right_pos];
                    }
                }
            }
        }

        // insert D3 and D5 positions
        vector<seq_len_t> D35_poses;
        D35_poses.reserve(seq_row_nonzeros + seq_col_nonzeros + 2);

        // TODO: start from zero?
        for (seq_len_t i = 0; i < seq_row_indices.size(); ++i) {
            for (seq_len_t j = 0; j < seq_row_num[seq_row_indices[i]]; ++j) {
                D35_poses.push_back(seq_row_indices[i]);
            }
        }
        for (seq_len_t i = 0; i < seq_col_indices.size(); ++i) { D35_poses.push_back(seq_col_indices[i]); }

        // Note! insert diversity gene seq poses BEFORE joining gene seq poses
        seq_poses.reserve(seq_poses.size() + D35_poses.size() + 2);  // +2 -> just in case (:
        seq_poses.insert(seq_poses.begin() + probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX), D35_poses.begin(), D35_poses.end());
    }


    void MAAGBuilder::buildNucVJinsertions(const ClonotypeNuc &clonotype, ProbMMC &probs, EventIndMMC &events,
                                           const vector<seq_len_t> &seq_poses, MetadataMode metadata_mode,
                                           ErrorMode error_mode) const {
        seq_len_t v_vertices = probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                j_vertices = probs.nodeRows(JOINING_DELETIONS_VJ_MATRIX_INDEX);

        probs.initNode(VarJoi_INSERTIONS_MATRIX_INDEX, 1, v_vertices, j_vertices);

        if (metadata_mode) {
            events.initNode(VarJoi_INSERTIONS_MATRIX_INDEX, 1, v_vertices, j_vertices);
        }

        this->buildNucInsertions(clonotype.sequence(),
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
                              *_vj_ins,
                              false);
    }


    void MAAGBuilder::buildNucVDinsertions(const ClonotypeNuc &clonotype, ProbMMC &probs, EventIndMMC &events,
                                           const vector<seq_len_t> &seq_poses, MetadataMode metadata_mode,
                                           ErrorMode error_mode) const {
        seq_len_t v_vertices = probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                d3_vertices = probs.nodeRows(DIVERSITY_GENES_MATRIX_INDEX);

        probs.initNode(VarDiv_INSERTIONS_MATRIX_INDEX, 1, v_vertices, d3_vertices);

        if (metadata_mode) {
            events.initNode(VarDiv_INSERTIONS_MATRIX_INDEX, 1, v_vertices, d3_vertices);
        }

        this->buildNucInsertions(clonotype.sequence(),
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
                              *_vd_ins,
                              false);
    }


    void MAAGBuilder::buildNucDJinsertions(const ClonotypeNuc &clonotype, ProbMMC &probs, EventIndMMC &events,
                                           const vector<seq_len_t> &seq_poses, MetadataMode metadata_mode,
                                           ErrorMode error_mode) const {
        seq_len_t v_vertices = probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX),
                d3_vertices = probs.nodeRows(DIVERSITY_GENES_MATRIX_INDEX),
                d5_vertices = probs.nodeColumns(DIVERSITY_GENES_MATRIX_INDEX),
                j_vertices = probs.nodeRows(JOINING_DELETIONS_VDJ_MATRIX_INDEX);

        probs.initNode(DivJoi_INSERTIONS_MATRIX_INDEX, 1, d5_vertices, j_vertices);

        if (metadata_mode) {
            events.initNode(DivJoi_INSERTIONS_MATRIX_INDEX, 1, d5_vertices, j_vertices);
        }

        this->buildNucInsertions(clonotype.sequence(),
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
                              *_dj_ins,
                              true);
    }


    void MAAGBuilder::buildNucInsertions(const string &sequence, ProbMMC &probs, EventIndMMC &events,
                                         const vector<seq_len_t> &seq_poses, ProbMMC::node_ind_t ins_node_index,
                                         event_ind_t null_insertion, seq_len_t max_size, MetadataMode metadata_mode,
                                         bool error_mode, seq_len_t left_vertices_start, seq_len_t left_vertices_end,
                                         seq_len_t right_vertices_start, seq_len_t right_vertices_end,
                                         const AbstractInsertionModel &mc, bool reversed) const {
        int insertion_len;
        bool good_insertion;
        char last_char = NULL_CHAR;

        for (size_t left_vertex_i = left_vertices_start; left_vertex_i <= left_vertices_end; ++left_vertex_i) {
            for (size_t right_vertex_i = right_vertices_start; right_vertex_i <= right_vertices_end; ++right_vertex_i) {

//                std::cout << (int) left_vertex_i << std::endl;
//                std::cout << (int) right_vertex_i << std::endl;
//                for (int i = 0; i <= seq_poses.size(); ++i) {
//                    std::cout << (int) seq_poses[i] << " ";
//                }
//                std::cout << std::endl;
//                for (int i = left_vertices_start; i <= left_vertices_end; ++i) {
//                    std::cout << (int) seq_poses[i] << " ";
//                }
//                std::cout << std::endl;
//                for (int i = right_vertices_start; i <= right_vertices_end; ++i) {
//                    std::cout << (int) seq_poses[i] << " ";
//                }
//                std::cout << std::endl;

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

//                std::cout << "done" << std::endl;
            }
        }
    }


    void MAAGBuilder::buildAAVariable(const ClonotypeAA &clonotype, ProbMMC &probs, CodonMMC &codons,
                                      vector<seq_len_t> &seq_poses) const
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
        codons.initNode(0, v_num, 1, len + 1);

        if (clonotype.recombination() == VJ_RECOMB) {
            probs.initNode(VARIABLE_GENES_MATRIX_INDEX, 1, v_num, j_num);
        } else if (clonotype.recombination() == VDJ_RECOMB) {
            probs.initNode(VARIABLE_GENES_MATRIX_INDEX, v_num, 1, 1);
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

            for (seq_len_t i = 1; i < v_end - v_start + 2; ++i) {
                codons(0, v_index, 0, i) = clonotype.getVarCodon(v_index, i);
            }
        }

        for (seq_len_t i = 0; i <= len; ++i) { seq_poses.push_back(i); }
    }


    void MAAGBuilder::buildAAJoining(const ClonotypeAA &clonotype, ProbMMC &probs, CodonMMC &codons,
                                     vector<seq_len_t> &seq_poses) const
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
        codons.initNode(codons.chainSize() - 1, j_num, len + 1, 1);

        // add J or J-D gene nodes
        if (clonotype.recombination() == VDJ_RECOMB) {
            probs.initNode(J_index_genes, 1, j_num, clonotype.nDiv());
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

            for (seq_len_t i = 0; i < clonotype.getJoiLen(j_index); ++i) {
                codons(codons.chainSize() - 1, j_index, i + shift, 0) = clonotype.getJoiCodon(j_index, i + 1);
            }
        }

        for (seq_len_t i = 3 * clonotype.sequence().size() - len + 1; i <= 3 * clonotype.sequence().size() + 1; ++i) {
            seq_poses.push_back(i);
        }
    }


    void MAAGBuilder::buildAADiversity(const ClonotypeAA &clonotype, ProbMMC &probs, CodonMMC &codons,
                                       vector<seq_len_t> &seq_poses) const
    {
        seq_len_t min_D_len;

        seq_len_t seq_arr_size = 3*clonotype.sequence().size() + 1; // TODO: +2 ??? Because 1-based with the range [0, size+1]
        // vector seq_row -> 0 means no such index in the matrix (row-wise), 1 otherwise
        std::vector<seq_len_t> seq_row(seq_arr_size, 0), seq_row_num = seq_row;

        // vector seq_col -> 0 means no such index in the matrix (column-wise), 1 otherwise.
        std::vector<seq_len_t> seq_col(seq_arr_size, 0);

        // shifts for each pair of positions (i,j) (in case if (i,j) is duplicated)
        std::vector<seq_len_t> shift_mat(seq_arr_size * seq_arr_size);

        // optimization for faster search for the max element in rows
        unique_ptr<bool[]> any_row_element(new bool[seq_arr_size]);

        for (seg_index_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
            min_D_len = _param_vec->D_min_len(clonotype.getDiv(d_index));
            std::fill(shift_mat.begin(), shift_mat.end(), 0);
            std::fill(any_row_element.get(), any_row_element.get() + seq_arr_size, false);

            for (seg_index_t d_align = 0; d_align < clonotype.numDivAlignments(d_index); ++d_align) {
                seq_len_t d_seq_start = clonotype.getDivSeqStart(d_index, d_align),
                        d_seq_end = clonotype.getDivSeqEnd(d_index, d_align);

                for (seq_len_t i = d_seq_start; i <= d_seq_end - min_D_len + 1; ++i) {
                    any_row_element[i] = true;

                    for (seq_len_t j = d_seq_start + min_D_len - static_cast<seq_len_t>(1); j <= d_seq_end; ++j) {
                        ++shift_mat[i*seq_arr_size + j];
                    }
                }

                for (seq_len_t j = d_seq_start + min_D_len - static_cast<seq_len_t>(1); j <= d_seq_end; ++j) {
                    seq_col[j] = 1;
                }
            }

            // for each row get max element, i.e., max number of repeats of a specific row index
            for (seq_len_t i = 0; i < seq_arr_size; ++i) {
                if (any_row_element[i]) {
                    seq_row[i] = std::max(seq_row[i],
                                          *std::max_element(shift_mat.begin() + i*seq_arr_size, shift_mat.begin() + (i+1)*seq_arr_size));
                }
            }
        }

        seq_len_t seq_row_nonzeros = 0, seq_ind = 0;
        std::vector<seq_len_t> seq_row_indices, seq_col_indices;
        seq_row_indices.reserve(seq_arr_size);
        seq_col_indices.reserve(seq_arr_size);

        // make new vector seq_start -> 1based index in rows of Ddel matrices
        for (seq_len_t i = 0; i < seq_arr_size; ++i) {
            if (seq_row[i]) {
                seq_row_indices.push_back(i);

                seq_row_num[i] = seq_row[i];
                seq_row_nonzeros += seq_row[i];

                seq_row[i] = seq_ind;
                seq_ind += seq_row_num[i];
            }
        }

        // make new vector seq_end -> 1based index in columns of Ddel matrices
        seq_ind = 0;
        for (seq_len_t i = 0; i < seq_arr_size; ++i) {
            if (seq_col[i]) {
                seq_col_indices.push_back(i);

                seq_col[i] = seq_ind;
                ++seq_ind;
            }
        }
        seq_len_t seq_col_nonzeros = seq_ind;


        probs.initNode(DIVERSITY_GENES_MATRIX_INDEX, clonotype.nDiv(), seq_row_nonzeros, seq_col_nonzeros);
        codons.initNode(1, clonotype.nDiv(), seq_row_nonzeros, seq_col_nonzeros); // start codons
        codons.initNode(2, clonotype.nDiv(), seq_row_nonzeros, seq_col_nonzeros); // end codons


        seg_index_t d_index, d_gene;
        seq_len_t d_len;
        seq_len_t d_seq_start, d_seq_end, d_gene_start, d_gene_end;
        codon_hash left_codon, right_codon;

        for (seg_index_t d_index = 0; d_index < clonotype.nDiv(); ++d_index) {
            d_gene = clonotype.getDiv(d_index);
            d_len = _genes->D()[d_gene].sequence.size();
            min_D_len = _param_vec->D_min_len(d_gene);

            std::fill(shift_mat.begin(), shift_mat.end(), 0);

            // for each aligned Div segment get all possible smaller alignments and add them to the matrix.
            for (seg_index_t j = 0; j < clonotype.numDivAlignments(d_index); ++j) {
                d_seq_start = clonotype.getDivSeqStart(d_index, j);
                d_seq_end = clonotype.getDivSeqEnd(d_index, j);
                d_gene_start = clonotype.getDivGeneStart(d_index, j);
                d_gene_end = clonotype.getDivGeneEnd(d_index, j);

                for (seq_len_t left_pos = d_seq_start; left_pos <= d_seq_end - min_D_len + 1; ++left_pos) {
//                    std::cout << (int) d_index << ":" << (int) j << std::endl;
//                    std::cout << d_seq_start << ":" << d_seq_end << ":" << d_gene_start << ":" << std::endl;
//                    std::cout << "left " << (int) left_pos << std::endl;

                    left_codon = clonotype.getDivCodon(d_index, j, left_pos - d_seq_start + 1);

                    for (seq_len_t right_pos = left_pos + min_D_len - 1; right_pos <= d_seq_end; ++right_pos) {
//                        std::cout << "right " << (int) right_pos << std::endl;
//                        std::cout << std::bitset<6>(left_codon).to_string() << std::endl;

                        right_codon = clonotype.getDivCodon(d_index, j, right_pos - left_pos + 1);

//                        std::cout << "---" << std::endl;
//                        std::cout << (int)left_pos << std::endl;
//                        std::cout << (int)right_pos << std::endl;
//                        std::cout << (int) (left_pos - d_seq_start + 1) << std::endl;
//                        std::cout << (int) (right_pos - left_pos + 1) << std::endl;
//                        std::cout << (int)left_codon << std::endl;
//                        std::cout << (int)right_codon << std::endl;
//                        std::cout << std::bitset<6>(right_codon).to_string() << std::endl;
//                        std::cout << "sizes " << (int) seq_row_nonzeros << ":" << (int) seq_col_nonzeros << std::endl;

                        probs(DIVERSITY_GENES_MATRIX_INDEX,
                              d_index,
                              seq_row[left_pos] + shift_mat[seq_arr_size*left_pos + right_pos],
                              seq_col[right_pos])
                                = _param_vec->event_prob(VDJ_DIV_DEL,
                                                         d_gene - 1,
                                                         d_gene_start + left_pos - d_seq_start,
                                                         d_len - (d_gene_end - (d_seq_end - right_pos)));

                        codons(1,
                               d_index,
                               seq_row[left_pos] + shift_mat[seq_arr_size*left_pos + right_pos],
                               seq_col[right_pos])
                                = left_codon;

                        codons(2,
                               d_index,
                               seq_row[left_pos] + shift_mat[seq_arr_size*left_pos + right_pos],
                               seq_col[right_pos])
                                = right_codon;

                        ++shift_mat[seq_arr_size*left_pos + right_pos];
                    }
                }
            }
        }

        // insert D3 and D5 positions
        vector<seq_len_t> D35_poses;
        D35_poses.reserve(seq_row_nonzeros + seq_col_nonzeros + 2);

        // TODO: start from zero?
        for (seq_len_t i = 0; i < seq_row_indices.size(); ++i) {
            for (seq_len_t j = 0; j < seq_row_num[seq_row_indices[i]]; ++j) {
                D35_poses.push_back(seq_row_indices[i]);
            }
        }
        for (seq_len_t i = 0; i < seq_col_indices.size(); ++i) { D35_poses.push_back(seq_col_indices[i]); }

        // Note! insert diversity gene seq poses BEFORE joining gene seq poses
        seq_poses.reserve(seq_poses.size() + D35_poses.size() + 2);  // +2 -> just in case (:
        seq_poses.insert(seq_poses.begin() + probs.nodeColumns(VARIABLE_DELETIONS_MATRIX_INDEX), D35_poses.begin(), D35_poses.end());
    }

}

#endif //YMIR_MAAG_BUILDER_H
