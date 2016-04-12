//
// Created by Vadim N. on 09/04/2016.
//

#ifndef YMIR_MAAG_BUILDER_BASE_H
#define YMIR_MAAG_BUILDER_BASE_H


#include "clonotype.h"
#include "genesegment.h"
#include "insertionmodel.h"
#include "maag.h"
#include "modelparametervector.h"
#include "repertoire.h"


namespace ymir {

    template <typename maag_type>
    class MAAGBuilderBase {

    public:

//        static const size_t verbose_step = 10000;


        /**
         * \brief Constructor for the builder from given vector with event probabilities and gene segments.
         */
        MAAGBuilderBase(const ModelParameterVector &param_vec, const VDJRecombinationGenes &genes)
                : MAAG(),
                  _param_vec(new ModelParameterVector(param_vec)),
                  _genes(new VDJRecombinationGenes(genes)) {
        }


        virtual ~MAAGBuilder()
        {
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

        unique_ptr<ModelParameterVector> _param_vec;
        unique_ptr<VDJRecombinationGenes> _genes;


        /**
         * \brief Private default constructor.
         */
        MAAGBuilderBase() : _param_vec(nullptr), _genes(nullptr) {}

}

#endif //YMIR_MAAG_BUILDER_BASE_H
