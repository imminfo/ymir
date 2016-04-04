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

#ifndef _ASSEMBLINGSTATISTICALMODEL_H_
#define _ASSEMBLINGSTATISTICALMODEL_H_


#include <unordered_set>

#include "maagbuilder.h"
#include "repertoire.h"
#include "clonotypeassembler.h"
#include "model_parser.h"


namespace ymir {

    /**
    * \class ProbabilisticAssemblingModel
    */
    class ProbabilisticAssemblingModel {

    public:

        /**
         * \brief Model constructor from the given folder with a model's parameters.
         *
         * \param folderpath Path to a folder with a model's parameters.
         * \param behav Model's behaviour. PREDEFINED means that each probability will stay constant for each input data
         * unless the model probabilities vector updated with new parameters (either all parameters or gene usage only). EMPTY means that the model will load
         * only gene segments sequences and generate vector of uniformly distributed marginal probabilities. This model could be used as an input model
         * for the statistical inference of the model's parameters.
         */
        ProbabilisticAssemblingModel(const string& folderpath, ModelBehaviour behav = PREDEFINED) {
            _model_path = folderpath + "/";
            _behaviour = behav;

            _recomb = UNDEF_RECOMB;
            _status = false;    

            this->parseModelConfig(_model_path);
            if (_status) {
                switch (_recomb) {
                    case VJ_RECOMB:
                        _parser.reset(new VJModelParser(_model_path, _config, behav));
                        break;
                    case VDJ_RECOMB:
                        _parser.reset(new VDJModelParser(_model_path, _config, behav));
                        break;
                    default: ;
                }

                _status = _parser->parse();
            }

            if (_status) {
                cout << "Model is loaded successfully." << endl;
                _parser->swap_genes(_genes);
                _parser->swap_parameters(_param_vec);
                this->make_builder();
                this->make_assembler();
            } else {
                cout << "Model is not loaded due to errors." << endl;
            }
        }


        /**
         *
         */
        ~ProbabilisticAssemblingModel()
        {
        }


        /**
         * \brief Compute full assembling probabilities of clonotypes from the given repertoire.
         * If some clonotypes have few alignments of some gene segments, than
         * full probability of all combinations of gene segments will be computed and either summed
         * or the max full probability will be chosen.
         *
         * \param repertoire Repertoire with alignments and event indices.
         * \param aminoacid What probabilities will be computed - of nucleotide or amino acid sequences. If
         * this parameter is false and some clone will have amino acid sequence, than error message will be generated to the output
         * and probability of this clone will set to zero.
         *
         * \return Vector of full assembling probabilities.
         */
        vector<prob_t> computeFullProbabilities(const ClonesetView& repertoire,
                                                ErrorMode error_mode,
                                                SequenceType sequence_type = NUCLEOTIDE,
                                                MAAGComputeProbAction action = SUM_PROBABILITY) const {
            return _builder->buildAndCompute(repertoire, error_mode, sequence_type, action);
        }


        /**
         * \brief Build a set of MAAGs from the given cloneset.
         *
         * \param repertoire Cloneset.
         * \param save_metadatae If true than build full MAAG (with additional matrices for event indices).
         * \param aminoacid If true than build graphs of aminoacid CDR3 sequences, traslating the given repertoire first.
         *
         * \return Set of MAAGs.
         */
        MAAGRepertoire buildGraphs(const ClonesetView &repertoire,
                                   MetadataMode save_metadata,
                                   ErrorMode error_mode,
                                   SequenceType sequence_type = NUCLEOTIDE,
                                   bool verbose = true) const {
#ifndef DNDEBUG
            if (!_status) {
                throw(std::runtime_error("Can't build graphs due to a model's failed status!"));
            }
#endif
            return _builder->build(repertoire, save_metadata, error_mode, sequence_type, verbose);
        }


        /**
         * \brief Update event probabilities in the given MAAG repertoire with new ones.
         */
        void updateEventProbabilities(MAAGRepertoire *repertoire, bool verbose = true) {
            _builder->updateEventProbabilities(repertoire, verbose);
        }


        /**
         * \brief Generate artificial repertoire of sequences using this model.
         *
         * \param count Number of generated sequences.
         * \param merge If true than return only unique clonotypes with counts of occurrences for each clonotype.
         *
         * \return Artificial repertoire.
         */
        Cloneset generateSequences(size_t count = 1) const {
#ifndef DNDEBUG
            if (!_status) {
                throw(std::runtime_error("Can't generate sequences due to a model's failed status!"));
            }
#endif
            // generate sequences
            return _generator->generate(count);
        }


        /**
         * \brief Access to the gene segments table.
         *
         * \return Const reference to the gene segments table.
         */
        const VDJRecombinationGenes& gene_segments() const {
#ifndef DNDEBUG
            if (!_status) {
                throw(std::runtime_error("Can't access gene segments in a model due to its failed status!"));
            }
#endif
            return *_genes;
        }


        /**
         * \brief Access to vector of probabilities.
         */
        const ModelParameterVector& event_probabilities() const {
#ifndef DNDEBUG
            if (!_status) {
                throw(std::runtime_error("Can't access the event probabilities vector in a model due to its failed status!"));
            }
#endif
            return *_param_vec;
        }

        /**
         * \brief Set a new event probabilities vector to this model.
         *
         * \param vec Vector with new event probabilities.
         */
        void updateModelParameterVector(const ModelParameterVector &vec) {
#ifndef DNDEBUG
            if (_param_vec->recombination() != vec.recombination()) {
                throw(std::runtime_error("Model parameter vectors are not comparable due to different recombination types!"));
            }

            if (_param_vec->size() != vec.size()) {
                throw(std::runtime_error("Model parameter vectors are not comparable due to different sizes!"));
            }
#endif

            *_param_vec = vec;
            _builder->updateModelParameterVector(vec);
            _generator->updateModelParameterVector(vec);
        }


        /**
         * \brief Given the cloneset, compute its gene usage on out-of-frames and
         * update the event probability vector with new gene probabilities.
         */
        void updateGeneUsage(const ClonesetView &cloneset) {
            prob_t laplace = 0;
            ClonesetView nonc = cloneset.noncoding();
            if (_recomb == VJ_RECOMB) {
                // Update V-J
                _param_vec->familyFill(VJ_VAR_JOI_GEN, 0);
                for (size_t i = 0; i < nonc.size(); ++i) {
                    uint n_vj = nonc[i].nVar() * nonc[i].nJoi();
                    for (seg_index_t v_i = 0; v_i < nonc[i].nVar(); ++v_i) {
                        for (seg_index_t j_i = 0; j_i < nonc[i].nJoi(); ++j_i) {
                            (*_param_vec)[_param_vec->event_index(VJ_VAR_JOI_GEN,
                                                                  0,
                                                                  nonc[i].getVar(v_i) - 1,
                                                                  nonc[i].getJoi(j_i) - 1)] += 1. / n_vj;
                        }
                    }
                }
                _param_vec->normaliseEventClass(VJ_VAR_JOI_GEN);
            } else if (_recomb == VDJ_RECOMB) {
                // Update V
                _param_vec->familyFill(VDJ_VAR_GEN, 0);
                for (size_t i = 0; i < nonc.size(); ++i) {
                    for (seg_index_t v_i = 0; v_i < nonc[i].nVar(); ++v_i) {
                        (*_param_vec)[_param_vec->event_index(VDJ_VAR_GEN,
                                                              0,
                                                              nonc[i].getVar(v_i) - 1)] += 1. / nonc[i].nVar();
                    }
                }
                _param_vec->normaliseEventClass(VDJ_VAR_GEN);

                // Update J-D
                _param_vec->familyFill(VDJ_JOI_DIV_GEN, 0);
                for (size_t i = 0; i < nonc.size(); ++i) {
                    uint n_jd = nonc[i].nJoi() * nonc[i].nDiv();
                    for (seg_index_t j_i = 0; j_i < nonc[i].nJoi(); ++j_i) {
                        for (seg_index_t d_i = 0; d_i < nonc[i].nDiv(); ++d_i) {
                            (*_param_vec)[_param_vec->event_index(VDJ_JOI_DIV_GEN,
                                                                  0,
                                                                  nonc[i].getJoi(j_i) - 1,
                                                                  nonc[i].getDiv(d_i) - 1)] += 1. / n_jd;
                        }
                    }
                }
                _param_vec->normaliseEventClass(VDJ_JOI_DIV_GEN);
            }
        }


        /**
        * \brief Save the model to the given folder.
        *
        * \param folderpath Path to the model folder.
        *
        * \return True if all is ok.
        */
        bool save(const string& folderpath) const {
            ofstream ofs;
            ofs.open(folderpath + "/model.json");

            if (ofs.is_open()) {

                ofs << _config;
                ofs.close();

                if (_recomb == VJ_RECOMB) {
                    this->save_vj(folderpath);
                } else if (_recomb == VDJ_RECOMB) {
                    this->save_vdj(folderpath);
                } else {
                    std::cout << "[ERROR] Can't save a model with an undefined recombination type." << std::endl;
                    return false;
                }

                return true;
            }

            std::cout << "[ERROR] Problem with saving a .json model file: probably there is no such directory: " << folderpath << std::endl;
            return false;
        }


        /**
         * \brief Return if the model is in the ready-to-work state.
         *
         * \return True or false.
         */
        bool status() const { return _status; }


        Recombination recombination() const { return _recomb; }


        string name() const { return _config.get("name", "Nameless model").asString(); }

    protected:

        bool _status, _verbose;
        ModelBehaviour _behaviour;

        Json::Value _config;
        Recombination _recomb;
        string _model_path;
        unique_ptr<ModelParser> _parser;

        unique_ptr<VDJRecombinationGenes> _genes;
        unique_ptr<ModelParameterVector> _param_vec;

        unique_ptr<MAAGBuilder> _builder;
        unique_ptr<ClonotypeAssembler> _generator;

        seq_len_t _min_D_len;


        /**
         * \brief Private default constructor.
         */
        ProbabilisticAssemblingModel() {}


        /**
         * \brief Parse JSON file with model parameters.
         */
        void parseModelConfig(const string& folderpath) {
            string jsonpath = folderpath + "/model.json";

            ifstream ifs;
            ifs.open(jsonpath);
            _status = false;
            if (ifs.is_open()) {
                ifs >> _config;
                if (_config.get("recombination", "undefined").asString() == "VJ") {
                    _recomb = VJ_RECOMB;
                    _status = true;
                } else if (_config.get("recombination", "undefined").asString() == "VDJ") {
                    _recomb = VDJ_RECOMB;
                    _status = true;
                } else {
                    _recomb = UNDEF_RECOMB;
                    _status = false;
                }

                cout << "Probabilistic assembling model:\n\t" <<
                        _config.get("name", "Nameless model").asString() <<
                        "\n\t(" <<
                        _config.get("comment", "").asString() <<
                        ")\n\t" <<
                        "Path : " << folderpath <<
                        "\n\t" <<
                        "Specimen : " << _config.get("specimen", "hobbit").asString() <<
                        "\n\t" <<
                        _config.get("recombination", "undefined").asString() <<
                        "-recombination  |  " <<
                        ((_config.get("errors", 0).asDouble() != 0) ? "Sequence error model" : "No sequence error model") <<
                        endl << "\tFiles:"<< endl;

                _status = true;
            } else {
                std::cout << "[ERROR] Probabilistic assembling model error:" << endl << "\t config .json file not found at [" << jsonpath << "]" << std::endl;
            }
        }


        /**
         * \brief Make MAAGBuilder with this model's parameters (event probabilities and gene segments).
         */
        void make_builder() {
            _builder.reset(new MAAGBuilder(*_param_vec, *_genes));
        }


        /**
         * \brief Make ClonotypeAssembler with this model's parameters (event probabilities and gene segments).
         */
        void make_assembler() {
            _generator.reset(new ClonotypeAssembler(*_param_vec, *_genes));
        }


        void save_vj(const string &folderpath) const {
            AbstractTDContainer* container;

            _genes->write(folderpath + _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("file", "vsegments.txt").asString(),
                          folderpath + _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("file", "jsegments.txt").asString());

            // V-J
            container = new TDMatrix(true);
            container->addMetadata(_genes->J().max());
            container->addColumnName("V / J");
            for (auto i = 1; i <= _genes->V().max(); ++i){
                container->addRowName(_genes->V()[i].allele);
            }
            for (auto i = 1; i <= _genes->J().max(); ++i){
                container->addColumnName(_genes->J()[i].allele);
            }
            container->addDataVector(_param_vec->get_iterator(1),
                                     _param_vec->get_iterator(_param_vec->eventClassSize(VJ_VAR_JOI_GEN) + 1));
            container->write(folderpath + _config.get("probtables", Json::Value()).get("v.j", Json::Value()).get("file", .0).asString());
            delete container;

            // V del
            container = new TDVectorList(true);

            for (auto i = 0; i <= _genes->V().maxLength(); ++i) {
                container->addRowName(std::to_string(i));
            }
            container->addColumnName("V deletions");

            for (auto i = 1; i <= _genes->V().max(); ++i) {
                container->addColumnName(_genes->V()[i].allele);
                container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VJ_VAR_DEL, i - 1, 0)),
                                         _param_vec->get_iterator(_param_vec->event_index(VJ_VAR_DEL, i - 1, 0) + _param_vec->eventFamilySize(VJ_VAR_DEL, i - 1)));
            }

            container->write(folderpath + _config.get("probtables", Json::Value()).get("v.del", Json::Value()).get("file", .0).asString());

            delete container;

            // J del
            container = new TDVectorList(true);
            for (auto i = 0; i <= _genes->J().maxLength(); ++i) {
                container->addRowName(std::to_string(i));
            }
            container->addColumnName("J deletions");
            for (auto i = 1; i <= _genes->J().max(); ++i){
                container->addColumnName(_genes->J()[i].allele);
                container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VJ_JOI_DEL, i - 1, 0)),
                                         _param_vec->get_iterator(_param_vec->event_index(VJ_JOI_DEL, i - 1, 0) + _param_vec->eventFamilySize(VJ_JOI_DEL, i - 1)));
            }
            container->write(folderpath + _config.get("probtables", Json::Value()).get("j.del", Json::Value()).get("file", .0).asString());
            delete container;

            // VJ ins
            container = new TDVectorList(true);
            container->addColumnName("VJ ins len");
            container->addColumnName("Probability");
            for (auto i = 0; i < _param_vec->eventClassSize(VJ_VAR_JOI_INS_LEN); ++i) {
                container->addRowName(std::to_string(i));
            }
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VJ_VAR_JOI_INS_LEN, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VJ_VAR_JOI_INS_LEN, 0, 0) + _param_vec->eventFamilySize(VJ_VAR_JOI_INS_LEN, 0)));
            container->write(folderpath + _config.get("probtables", Json::Value()).get("ins.len", Json::Value()).get("file", .0).asString());
            delete container;

            // VJ nuc
            container = new TDVectorList(true);
            container->addColumnName("VJ nucs");
            container->addColumnName("Probability");
            container->addRowName("A"); container->addRowName("C"); container->addRowName("G"); container->addRowName("T");
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VJ_VAR_JOI_INS_NUC, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VJ_VAR_JOI_INS_NUC, 0, 0) + 4));
            container->write(folderpath + _config.get("probtables", Json::Value()).get("ins.nucl", Json::Value()).get("file", .0).asString());

            delete container;
        }


        void save_vdj(const string &folderpath) const {
            AbstractTDContainer* container;

            _genes->write(folderpath + _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("file", "vsegments.txt").asString(),
                          folderpath + _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("file", "jsegments.txt").asString(),
                          folderpath + _config.get("segments", Json::Value("")).get("diversity", Json::Value("")).get("file", "dsegments.txt").asString());

            // V
            container = new TDVector(true);
            container->addColumnName("V");
            container->addColumnName("Probability");
            container->addDataVector(vector<prob_t>());
            for (auto i = 1; i <= _genes->V().max(); ++i) {
                container->addRowName(_genes->V()[i].allele);
                container->addDataValue(_param_vec->event_prob(VDJ_VAR_GEN, 0, i - 1));
            }
            container->write(folderpath + _config.get("probtables", Json::Value()).get("v", Json::Value()).get("file", .0).asString());
            delete container;

            // J-D
            container = new TDMatrix(true);
            container->addMetadata(_genes->D().max());
            container->addColumnName("D / J");
            for (auto i = 1; i <= _genes->J().max(); ++i) {
                container->addRowName(_genes->J()[i].allele);
            }
            for (auto i = 1; i <= _genes->D().max(); ++i) {
                container->addColumnName(_genes->D()[i].allele);
            }
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_JOI_DIV_GEN, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->eventClassSize(VDJ_JOI_DIV_GEN) + _param_vec->event_index(VDJ_JOI_DIV_GEN, 0, 0)));
            container->write(folderpath + _config.get("probtables", Json::Value()).get("j.d", Json::Value()).get("file", .0).asString());
            delete container;

            // V del
            container = new TDVectorList(true);
            for (auto i = 0; i <= _genes->V().maxLength(); ++i) {
                container->addRowName(std::to_string(i));
            }
            container->addColumnName("V deletions");
            for (auto i = 1; i <= _genes->V().max(); ++i) {
                container->addColumnName(_genes->V()[i].allele);
                container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DEL, i - 1, 0)),
                                         _param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DEL, i - 1, 0) + _param_vec->eventFamilySize(VDJ_VAR_DEL, i - 1)));
            }
            container->write(folderpath + _config.get("probtables", Json::Value()).get("v.del", Json::Value()).get("file", .0).asString());
            delete container;

            // J del
            container = new TDVectorList(true);
            for (auto i = 0; i <= _genes->J().maxLength(); ++i) {
                container->addRowName(std::to_string(i));
            }
            container->addColumnName("J deletions");
            for (auto i = 1; i <= _genes->J().max(); ++i){
                container->addColumnName(_genes->J()[i].allele);
                container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_JOI_DEL, i - 1, 0)),
                                         _param_vec->get_iterator(_param_vec->event_index(VDJ_JOI_DEL, i - 1, 0) + _param_vec->eventFamilySize(VDJ_JOI_DEL, i - 1)));
            }
            container->write(folderpath + _config.get("probtables", Json::Value()).get("j.del", Json::Value()).get("file", .0).asString());
            delete container;

            // D del
            container = new TDMatrixList(true);
            for (auto i = 1; i <= _genes->D().max(); ++i) {
                container->addRowName(_genes->D()[i].allele);
                container->addMetadata(_param_vec->n_columns(VDJ_DIV_DEL, i - 1));
                container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_DEL, i - 1, 0, 0)),
                                         _param_vec->get_iterator(_param_vec->eventFamilySize(VDJ_DIV_DEL, i - 1) + _param_vec->event_index(VDJ_DIV_DEL, i - 1, 0, 0)));
            }
            container->write(folderpath + _config.get("probtables", Json::Value()).get("d.del", Json::Value()).get("file", .0).asString());
            delete container;

            // VD ins + DJ ins
            container = new TDVectorList(true);
            container->addColumnName("Insertion length");
            container->addColumnName("VD ins");
            container->addColumnName("DJ ins");
            for (auto i = 0; i < _param_vec->eventClassSize(VDJ_VAR_DIV_INS_LEN); ++i) {
                container->addRowName(std::to_string(i));
            }
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_LEN, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_LEN, 0, 0) + _param_vec->eventFamilySize(VDJ_VAR_DIV_INS_LEN, 0)));
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_LEN, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_LEN, 0, 0) + _param_vec->eventFamilySize(VDJ_DIV_JOI_INS_LEN, 0)));
            container->write(folderpath + _config.get("probtables", Json::Value()).get("ins.len", Json::Value()).get("file", .0).asString());
            delete container;

            // VD nuc + DJ nuc
            container = new TDVectorList(true);
            container->addColumnName("Prev / Next");
            container->addColumnName("VD.A"); container->addColumnName("VD.C"); container->addColumnName("VD.G"); container->addColumnName("VD.T");
            container->addColumnName("DJ.A"); container->addColumnName("DJ.C"); container->addColumnName("DJ.G"); container->addColumnName("DJ.T");
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC_A, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC_A, 0, 0) + 4));
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC_C, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC_C, 0, 0) + 4));
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC_G, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC_G, 0, 0) + 4));
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC_T, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_VAR_DIV_INS_NUC_T, 0, 0) + 4));
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC_A, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC_A, 0, 0) + 4));
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC_C, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC_C, 0, 0) + 4));
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC_G, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC_G, 0, 0) + 4));
            container->addDataVector(_param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC_T, 0, 0)),
                                     _param_vec->get_iterator(_param_vec->event_index(VDJ_DIV_JOI_INS_NUC_T, 0, 0) + 4));
            container->write(folderpath + _config.get("probtables", Json::Value()).get("ins.nucl", Json::Value()).get("file", .0).asString());
            delete container;
        }
    };

}


#endif