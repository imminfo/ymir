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

#include "textdata.h"
#include "maagbuilder.h"
#include "repertoire.h"
#include "genesegment.h"
#include "clonotypeassembler.h"


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

            _recomb = UNDEFINED;
            _status = false;

            _genes = nullptr;
            _param_vec = nullptr;

            _builder = nullptr;
            _generator = nullptr;

            _status = this->parseModelConfig(_model_path + "model.json")
                      && this->parseGeneSegments()
                      && this->makeEventProbabilitiesVector();

            if (_status) {
                cout << "Model is loaded successfully." << endl;
                this->make_builder();
                this->make_assembler();
            } else {
                cout << "Model is not loaded due to errors." << endl;
            }
        }


        /**
         *
         */
        virtual ~ProbabilisticAssemblingModel() {
            if (_genes) { delete _genes; }
            if (_param_vec) { delete _param_vec; }
            if (_builder) { delete _builder; }
            if (_generator) { delete _generator; }
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
                                                bool aminoacid = false,
                                                MAAGComputeProbAction action = SUM_PROBABILITY) const {
            return this->_builder->buildAndCompute(repertoire, aminoacid, action);
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
                                   MetadataMode save_metadata = SAVE_METADATA,
                                   SequenceType sequence_type = NUCLEOTIDE,
                                   bool verbose = true) const {
            return _builder->build(repertoire, save_metadata, verbose);
        }


        /**
         * \brief Update event probabilities in the given MAAG repertoire with new ones.
         */
        void updateEventProbabilities(MAAGRepertoire *repertoire) {
            this->_builder->updateEventProbabilities(repertoire);
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
            // generate sequences
            return _generator->generate(count);
        }


        /**
         * \brief Access to the gene segments table.
         *
         * \return Const reference to the gene segments table.
         */
        const VDJRecombinationGenes& gene_segments() const {
#ifdef YDEBUG
            if (_status) {
                throw(std::runtime_error("Can't access gene segments in a model due to its failed status!"));
            }
#endif
            return *_genes;
        }


        /**
         * \brief Access to vector of probabilities.
         */
        const ModelParameterVector& event_probabilities() const {
#ifdef YDEBUG
            if (_status) {
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
        void updateEventProbabilitiesVector(const ModelParameterVector &vec) {
#ifdef YDEBUG
            if (_param_vec->recombination() != vec.recombination()) {
                throw(std::runtime_error("Model parameter vectors are not comparable due to different recombination types!"));
            }

            if (_param_vec->size() != vec.size()) {
                throw(std::runtime_error("Model parameter vectors are not comparable due to different sizes!"));
            }
#endif

            *_param_vec = vec;
            this->make_builder();
            this->make_assembler();
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
                _param_vec->normaliseEventFamily(VJ_VAR_JOI_GEN);

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
                _param_vec->normaliseEventFamily(VDJ_VAR_GEN);

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
                _param_vec->normaliseEventFamily(VDJ_JOI_DIV_GEN);
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


        string name() const { return _config.get("name", "Nameless model").asString(); }

    protected:

        bool _status, _verbose;
        ModelBehaviour _behaviour;

        Json::Value _config;
        Recombination _recomb;
        string _model_path;

        VDJRecombinationGenes *_genes;
        ModelParameterVector *_param_vec;

        MAAGBuilder *_builder;
        ClonotypeAssembler *_generator;

        seq_len_t _min_D_len;


        /**
         * \brief Private default constructor.
         */
        ProbabilisticAssemblingModel() {}


        /**
         * \brief Parse JSON file with model parameters.
         */
        bool parseModelConfig(const string& jsonpath) {
            ifstream ifs;
            ifs.open(jsonpath);
            if (ifs.is_open()) {
                ifs >> _config;
                if (_config.get("recombination", "undefined").asString() == "VJ") {
                    _recomb = VJ_RECOMB;
                } else if (_config.get("recombination", "undefined").asString() == "VDJ") {
                    _recomb = VDJ_RECOMB;
                } else {
                    _recomb = UNDEFINED;
                }
                cout << "Probabilistic assembling model:\n\t" <<
                        _config.get("name", "Nameless model").asString() <<
                        "\n\t(" <<
                        _config.get("comment", "").asString() <<
                        ")\n\t" <<
                        _config.get("recombination", "undefined").asString() <<
                        "-recombination  |  " <<
                        (_config.get("hypermutations", false).asBool() ? "Hypermutations" : "No hypermutations") <<
                        endl << "\tFiles:"<< endl;
                return true;
            }
            std::cout << "[ERROR] Probabilistic assembling model error:" << endl << "\t config .json file not found at [" << jsonpath << "]" << std::endl;
            return false;
        }


        /**
         * \brief Parse gene segment JSON files and tables.
         */
        bool parseGeneSegments() {
            cout << "\tV gene seg.:     ";
            if (_config.get("segments", Json::Value("")).get("variable", Json::Value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;
                return false;
            } else { cout << "OK" << endl; }
            string v_path = _model_path + _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("file", "").asString();

            cout << "\tJ gene seg.:     ";
            if (_config.get("segments", Json::Value("")).get("joining", Json::Value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;
                return false;
            } else { cout << "OK" << endl; }
            string j_path = _model_path + _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("file", "").asString();

            bool vok, jok, dok = true;

            if (_recomb == VJ_RECOMB) {
                _genes = new VDJRecombinationGenes("VJ.V", v_path, "VJ.J", j_path, &vok, &jok);
                /* add P nucleotides */
            } else if (_recomb == VDJ_RECOMB) {
                cout << "\tD gene seg.:     ";
                if (_config.get("segments", Json::Value("")).get("diversity", Json::Value("")).size() == 0) {
                    cout << "ERROR: no gene segments file in the model's .json." << endl;
                    return false;
                } else {
                    string d_path = _model_path + _config.get("segments", Json::Value("")).get("diversity", Json::Value("")).get("file", "").asString();
                    _genes = new VDJRecombinationGenes("VDJ.V", v_path, "VDJ.J", j_path, "VDJ.D", d_path, &vok, &jok, &dok);
                    _min_D_len = _config.get("segments", Json::Value("")).get("diversity", Json::Value("")).get("min.len", DEFAULT_DIV_GENE_MIN_LEN).asUInt();
                    /* add P nucleotides */
                    cout << "OK" << endl;
                }
            } else {
                cout << "Undefined recombination type." << endl;
            }
            return vok && jok && dok;
        }


        /**
         * \brief Parse files with event probabilities matrices and make ModelParameterVector.
         */
        bool makeEventProbabilitiesVector() {
            if (_behaviour == EMPTY) {
                return this->createEventProbabilitiesFromScratch();
            } else {
                return this->parseEventProbabilitiesFromFiles();
            }
        }


        /**
         *
         */
        bool findGenes(const vector<string> &names, const GeneSegmentAlphabet &gsa, string &err_message) const {
            unordered_set<string> nameset;

            for (size_t i = 0; i < names.size(); ++i) {
                if (gsa[names[i]].index == 0) {
                    err_message = "ERROR: can't find " + names[i] + " in gene segments.";
                    return false;
                }
                nameset.insert(names[i]);
            }

            for (size_t i = 1; i < gsa.size(); ++i) {
                if (nameset.count(gsa[i].allele) == 0) {
                    err_message = "ERROR: can't find " + gsa[i].allele + " in this file.";
                    return false;
                }
            }

            return true;
        }


        /**
         *
         */
        vector<seg_index_t> arrangeNames(const vector<string> &names, const GeneSegmentAlphabet &gsa) const {
            vector<seg_index_t> res;
            res.resize(names.size(), 0);

            for (size_t i = 0; i < names.size(); ++i) { res[i] = gsa[names[i]].index; }

            return res;
        }


        void addGenes(AbstractTDContainer *container,
                      const GeneSegmentAlphabet &gsa,
                      vector<prob_t> &event_probs,
                      vector<event_ind_t> &event_lengths,
                      vector<event_ind_t> &event_classes,
                      vector<seq_len_t> &event_col_num,
                      vector<prob_t> &laplace) const {
            vector<seg_index_t> name_order = this->arrangeNames(container->row_names(), gsa);
            vector<prob_t> prob_data;
            prob_data.resize(container->data(0).size(), 0);
            for (size_t i = 0; i < name_order.size(); ++i) {
                prob_data[name_order[i] - 1] = container->data(0)[i];
            }
            event_probs.insert(event_probs.end(),
                               prob_data.begin(),
                               prob_data.end());
            event_lengths.push_back(prob_data.size());
            event_col_num.push_back(0);
            laplace.push_back(container->laplace());

            event_classes.push_back(0);
        }


        void addGenes(AbstractTDContainer *container,
                      const GeneSegmentAlphabet &gsa_row,
                      const GeneSegmentAlphabet &gsa_column,
                      vector<prob_t> &event_probs,
                      vector<event_ind_t> &event_lengths,
                      vector<event_ind_t> &event_classes,
                      vector<seq_len_t> &event_col_num,
                      vector<prob_t> &laplace,
                      seg_index_t prev_class_size) const {
            vector<seg_index_t> name_order_row = this->arrangeNames(container->row_names(), gsa_row);
            vector<seg_index_t> name_order_column = this->arrangeNames(container->column_names(), gsa_column);
            vector<prob_t> prob_data = container->data(0), sorted_prob_data;
            sorted_prob_data.resize(prob_data.size(), 0);
            for (size_t i = 0; i < name_order_row.size(); ++i) {
                for (size_t j = 0; j < name_order_column.size(); ++j) {
                    sorted_prob_data[(name_order_row[i] - 1) * container->n_columns() + (name_order_column[j] - 1)] = prob_data[i * container->n_columns() + j];
                }
            }
            event_probs.insert(event_probs.end(),
                               sorted_prob_data.begin(),
                               sorted_prob_data.end());
            event_lengths.push_back(sorted_prob_data.size());
            event_col_num.push_back(container->n_columns());
            laplace.push_back(container->laplace());

            if (prev_class_size) {
                event_classes.push_back(event_classes[event_classes.size() - 1] + prev_class_size);
            } else {
                event_classes.push_back(0);
            }
        }


        void addDels(AbstractTDContainer *container,
                     const GeneSegmentAlphabet &gsa,
                     vector<prob_t> &event_probs,
                     vector<event_ind_t> &event_lengths,
                     vector<event_ind_t> &event_classes,
                     vector<seq_len_t> &event_col_num,
                     vector<prob_t> &laplace,
                     seg_index_t prev_class_size) const {
            vector<seg_index_t> name_order = this->arrangeNames(container->column_names(), gsa);
            vector<prob_t> prob_data;
            for (size_t i = 0; i < name_order.size(); ++i) {
                // find correct segment for i-th position
                size_t j = 0;
                for (; i+1 != name_order[j] ; ++j) {}
                prob_data = container->data(j);

                // add trailing zeros if distribution is smaller than a gene length
                if (prob_data.size() < gsa[name_order[j]].sequence.size() + 1) {
                    prob_data.resize(gsa[name_order[j]].sequence.size() + 1, 0);
                }

                // remove trailing zeros
                if (prob_data.size() > gsa[name_order[j]].sequence.size() + 1) {
                    prob_data.resize(gsa[name_order[j]].sequence.size() + 1);
                }

                event_probs.insert(event_probs.end(),
                                   prob_data.begin(),
                                   prob_data.end());
                event_lengths.push_back(prob_data.size());
                event_col_num.push_back(0);
                laplace.push_back(container->laplace());
            }
            event_classes.push_back(event_classes[event_classes.size() - 1] + prev_class_size);
        }


        void addDels2D(AbstractTDContainer *container,
                       const GeneSegmentAlphabet &gsa,
                       vector<prob_t> &event_probs,
                       vector<event_ind_t> &event_lengths,
                       vector<event_ind_t> &event_classes,
                       vector<seq_len_t> &event_col_num,
                       vector<prob_t> &laplace,
                       seg_index_t prev_class_size) const {
            vector<seg_index_t> name_order = this->arrangeNames(container->row_names(), gsa);
            vector<prob_t> prob_data;
            for (size_t i = 0; i < name_order.size(); ++i) {
                // find correct segment for i-th position
                size_t j = 0;
                for (; i+1 != name_order[j] ; ++j) {}
                prob_data = container->data(j);

                if (prob_data.size() > gsa[name_order[j]].sequence.size() + 1) {
                    vector<prob_t> new_prob_data((gsa[name_order[j]].sequence.size() + 1) * (gsa[name_order[j]].sequence.size() + 1));
                    for (auto row_i = 0; row_i < gsa[name_order[j]].sequence.size() + 1; ++row_i) {
                        for (auto col_i = 0; col_i < gsa[name_order[j]].sequence.size() + 1; ++col_i) {
                            new_prob_data[row_i * (gsa[name_order[j]].sequence.size() + 1) + col_i] = prob_data[row_i * (gsa[name_order[j]].sequence.size() + 1) + col_i];
                        }
                    }
                    prob_data = new_prob_data;
                }

                event_probs.insert(event_probs.end(),
                                   prob_data.begin(),
                                   prob_data.end());
                event_lengths.push_back(prob_data.size());
                event_col_num.push_back(container->metadata(j));
                laplace.push_back(container->laplace());
            }
            event_classes.push_back(event_classes[event_classes.size() - 1] + prev_class_size);
        }


        void addIns(AbstractTDContainer *container,
                    vector<prob_t> &event_probs,
                    vector<event_ind_t> &event_lengths,
                    vector<event_ind_t> &event_classes,
                    vector<seq_len_t> &event_col_num,
                    vector<prob_t> &laplace,
                    seg_index_t prev_class_size,
                    seq_len_t max_ins_len = 0) const {
            vector<prob_t> prob_data;
            for (size_t i = 0; i < container->n_columns(); ++i) {
                prob_data = container->data(i);
                if (max_ins_len) { prob_data.resize(max_ins_len); }
                event_probs.insert(event_probs.end(),
                                   prob_data.begin(),
                                   prob_data.end());
                event_lengths.push_back(prob_data.size());
                event_col_num.push_back(0);
                laplace.push_back(container->laplace());
                event_classes.push_back(event_classes[event_classes.size() - 1] + prev_class_size);
                prev_class_size = 1;
            }
        }


        /**
         * \brief Make MAAGBuilder with this model's parameters (event probabilities and gene segments).
         */
        void make_builder() {
            if (_builder) { delete _builder; }
            _builder = new MAAGBuilder(*_param_vec, *_genes);
        }


        /**
         * \brief Make ClonotypeAssembler with this model's parameters (event probabilities and gene segments).
         */
        void make_assembler() {
            if (_generator) { delete _generator; }
            _generator = new ClonotypeAssembler(*_param_vec, *_genes);
        }


        /**
         *
         */
        bool parseEventProbabilitiesFromFiles() {
            Json::Value pt = _config.get("probtables", "no-prob");
            if (pt.size()) {
                AbstractTDContainer *container = nullptr;
                string element = "", err_message = "";
                vector<AbstractTDContainer*> containers;
                containers.resize(10, nullptr);

                // Parser tables with event probabilities
                for (Json::ArrayIndex i = 0; i < pt.size(); ++i) {
                    element = pt.getMemberNames()[i];
                    container = read_textdata(_model_path + pt[element]["file"].asString(),
                                              pt[element]["type"].asString(),
                                              pt[element].get("skip.first.column", true).asBool(),
                                              pt[element].get("laplace", .0).asDouble(),
                                              err_message);

                    this->parseDataContainer(element, container, containers);
                }

                return this->makeModelParameterVector(containers);
            }

            std::cout << "[ERROR] No information about probability events in the model .json file found." << std::endl;
            return false;
        }


        void parseDataContainer(const string &element, AbstractTDContainer *container, vector<AbstractTDContainer*> &containers) {
            if (_recomb == VJ_RECOMB) {
                this->parseVJContainer(element, container, containers);
            }
            else if (_recomb == VDJ_RECOMB) {
                this->parseVDJContainer(element, container, containers);
            }
        }


        void parseVJContainer(const string &element, AbstractTDContainer *container, vector<AbstractTDContainer*> &containers) {
            string err_message = "OK";
            if (element == "v.j") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->J(), err_message)
                    && this->findGenes(container->row_names(), _genes->V(), err_message)) {
                    containers[VJ_VAR_JOI_GEN] = container;
                }

                cout << "\tV-J gene pairs:  " << err_message << endl;
            }
            else if (element == "v.del") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->V(), err_message)) {
                    containers[VJ_VAR_DEL] = container;
                }

                cout << "\tV delet. num.:   " << err_message << endl;
            }
            else if (element == "j.del") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->J(), err_message)) {
                    containers[VJ_JOI_DEL] = container;
                }

                cout << "\tJ delet. num.:   " << err_message << endl;
            }
            else if (element == "ins.len") {
                if (container && container->file_exists()) {
                    if (container->n_columns() != 1) {
                        stringstream ss;
                        ss << "ERROR: wrong number of columns (expected: 1, got: " << (int) container->n_columns() << ")";
                        err_message = ss.str();
                    } else {
                        containers[VJ_VAR_JOI_INS_LEN] = container;
                    }
                }

                cout << "\tVJ ins. len.:    " << err_message << endl;
            }
            else if (element == "ins.nucl") {
                if (container && container->file_exists()) {
                    if (container->n_rows() != 4 || container->n_columns() != 1) {
                        stringstream ss;
                        ss << "ERROR: wrong number of columns and rows (expected: 4 X 1, got: " << (int) container->n_rows() << "X" << (int) container->n_columns() << ")";
                        err_message = ss.str();
                    } else {
                        containers[VJ_VAR_JOI_INS_NUC] = container;
                    }
                }

                cout << "\tVJ ins. nuc.:    " << err_message << endl;
            }
            else { cerr << "Unrecognised element in \'probtables\'" << ":\n\t" << element << endl; }
        }


        void parseVDJContainer(const string &element, AbstractTDContainer *container, vector<AbstractTDContainer*> &containers) {
            string err_message = "OK";
            if (element == "v") {
                if (container && container->file_exists()) { containers[VDJ_VAR_GEN] = container; }
                cout << "\tV genes prob.:   " << err_message << endl;
            }
            else if (element == "j.d") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->D(), err_message)
                    && this->findGenes(container->row_names(), _genes->J(), err_message)) {
                    containers[VDJ_JOI_DIV_GEN] = container;
                }

                cout << "\tJ-D gene pairs:  " << err_message << endl;
            }
            else if (element == "v.del") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->V(), err_message)) {
                    containers[VDJ_VAR_DEL] = container;
                }

                cout << "\tV delet. num.:   " << err_message << endl;;
            }
            else if (element == "j.del") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->J(), err_message)) {
                    containers[VDJ_JOI_DEL] = container;
                }

                cout << "\tJ delet. num.:   " << err_message << endl;
            }
            else if (element == "d.del") {
                if (container && container->file_exists()) { containers[VDJ_DIV_DEL] = container; }

                cout << "\tD delet. num.:   " << err_message << endl;
            }
            else if (element == "ins.len") {
                if (container && container->file_exists()) {
                    if (container->n_columns() != 2) {
                        stringstream ss;
                        ss << "ERROR: wrong number of columns (expected: 2, got: " << (int) container->n_columns() << ")";
                        err_message = ss.str();
                    } else {
                        containers[VDJ_VAR_DIV_INS_LEN] = container;
                    }
                }

                cout << "\tVD/DJ ins. len.: " << err_message << endl;;
            }
            else if (element == "ins.nucl") {
                if (container && container->file_exists()) {
                    if (container->n_rows() != 4 || container->n_columns() != 8) {
                        stringstream ss;
                        ss << "ERROR: wrong number of columns and rows (expected: 4 X 8, got: " << (int) container->n_rows() << "X" << (int) container->n_columns() << ")";
                        err_message = ss.str();
                    } else {
                        containers[VDJ_VAR_DIV_INS_NUC] = container;
                    }
                }

                cout << "\tVD/DJ ins. nuc.: " << err_message << endl;
            }
            else { cerr << "Unrecognised element in \'probtables\'" << ":\n\t" << element << endl; }
        }


        /**
         * \brief Create marginal probabilities based on gene segments sequences and make all event probabilities uniform.
         */
        bool createEventProbabilitiesFromScratch() {
            vector<AbstractTDContainer*> containers;
            containers.resize(10, nullptr);

            if (_recomb == VJ_RECOMB) {
                this->createVJContainers(containers);
            } else if (_recomb == VDJ_RECOMB) {
                this->createVDJContainers(containers);
            } else {
                return false;
            }

            bool is_ok = this->makeModelParameterVector(containers);

            _param_vec->fill(1);
            _param_vec->normaliseEventFamilies();

            if (!is_ok) { throw(runtime_error("WRONG EMPTY VECTOR CREATING SUBROUTINE!!!")); }

            return is_ok;
        }


        /**
         * \brief Create new marginal probabilities for a VJ recombination generation model.
         */
        void createVJContainers(vector<AbstractTDContainer*> &containers) {
            AbstractTDContainer* container;

            // V-J
            container = new TDMatrix(true, _config.get("probtables", Json::Value()).get("v.j", Json::Value()).get("laplace", .0).asDouble());
            for (auto i = 1; i <= _genes->V().max(); ++i){
                container->addRowName(_genes->V()[i].allele);
            }
            for (auto i = 1; i <= _genes->J().max(); ++i){
                container->addColumnName(_genes->J()[i].allele);
            }
            container->addDataVector(vector<prob_t>(_genes->V().max() * _genes->J().max()));
            containers[VJ_VAR_JOI_GEN] = container;

            // V del
            container = new TDVectorList(true, _config.get("probtables", Json::Value()).get("v.del", Json::Value()).get("laplace", .0).asDouble());
            for (auto i = 1; i <= _genes->V().max(); ++i) {
                container->addColumnName(_genes->V()[i].allele);
                container->addDataVector(vector<prob_t>(_genes->V()[i].sequence.size() + 1));
            }
            containers[VJ_VAR_DEL] = container;

            // J del
            container = new TDVectorList(true, _config.get("probtables", Json::Value()).get("j.del", Json::Value()).get("laplace", .0).asDouble());
            for (auto i = 1; i <= _genes->J().max(); ++i){
                container->addColumnName(_genes->J()[i].allele);
                container->addDataVector(vector<prob_t>(_genes->J()[i].sequence.size() + 1));
            }
            containers[VJ_JOI_DEL] = container;

            // VJ ins
            container = new TDVectorList(true, _config.get("probtables", Json::Value()).get("ins.len", Json::Value()).get("laplace", .0).asDouble());
            container->addColumnName("VJ ins len");
            container->addDataVector(vector<prob_t>(_config.get("probtables", Json::Value()).get("ins.len", Json::Value()).get("max.len", DEFAULT_MAX_INS_LENGTH).asUInt64() + 1));
            containers[VJ_VAR_JOI_INS_LEN] = container;

            // VJ nuc
            container = new TDVectorList(true, _config.get("probtables", Json::Value()).get("ins.nucl", Json::Value()).get("laplace", .0).asDouble());
            container->addColumnName("VJ nucs");
            container->addDataVector(vector<prob_t>(4));
            containers[VJ_VAR_JOI_INS_NUC] = container;
        }


        /**
         * \brief Create new marginal probabilities for a VDJ recombination generation model.
         */
        void createVDJContainers(vector<AbstractTDContainer*> &containers) {
            AbstractTDContainer* container;

            // V
            container = new TDVector(true, _config.get("probtables", Json::Value()).get("v", Json::Value()).get("laplace", .0).asDouble());
            container->addDataVector(vector<prob_t>());
            for (auto i = 1; i <= _genes->V().max(); ++i) {
                container->addRowName(_genes->V()[i].allele);
                container->addDataValue(1);
            }
            containers[VDJ_VAR_GEN] = container;

            // J-D
            container = new TDMatrix(true, _config.get("probtables", Json::Value()).get("j.d", Json::Value()).get("laplace", .0).asDouble());
            for (auto i = 1; i <= _genes->J().max(); ++i) {
                container->addRowName(_genes->J()[i].allele);
            }
            for (auto i = 1; i <= _genes->D().max(); ++i) {
                container->addColumnName(_genes->D()[i].allele);
            }
            container->addDataVector(vector<prob_t>(_genes->D().max() * _genes->J().max()));
            containers[VDJ_JOI_DIV_GEN] = container;

            // V del
            container = new TDVectorList(true, _config.get("probtables", Json::Value()).get("v.del", Json::Value()).get("laplace", .0).asDouble());
            for (auto i = 1; i <= _genes->V().max(); ++i) {
                container->addColumnName(_genes->V()[i].allele);
                container->addDataVector(vector<prob_t>(_genes->V()[i].sequence.size() + 1));
            }
            containers[VDJ_VAR_DEL] = container;

            // J del
            container = new TDVectorList(true, _config.get("probtables", Json::Value()).get("j.del", Json::Value()).get("laplace", .0).asDouble());
            for (auto i = 1; i <= _genes->J().max(); ++i) {
                container->addColumnName(_genes->J()[i].allele);
                container->addDataVector(vector<prob_t>(_genes->J()[i].sequence.size() + 1));
            }
            containers[VDJ_JOI_DEL] = container;

            // D del
            container = new TDMatrixList(true, _config.get("probtables", Json::Value()).get("d.del", Json::Value()).get("laplace", .0).asDouble());
            for (auto i = 1; i <= _genes->D().max(); ++i) {
                container->addColumnName(_genes->D()[i].allele);
                container->addDataVector(vector<prob_t>( (_genes->D()[i].sequence.size() + 1) * (_genes->D()[i].sequence.size() + 1) ));
                container->addRowName(_genes->D()[i].allele);
                container->addMetadata(_genes->D()[i].sequence.size() + 1);
            }
            containers[VDJ_DIV_DEL] = container;

            // VD ins + DJ ins
            container = new TDVectorList(true, _config.get("probtables", Json::Value()).get("ins.len", Json::Value()).get("laplace", .0).asDouble());
            container->addColumnName("VD ins");
            container->addColumnName("DJ ins");
            container->addDataVector(vector<prob_t>(_config.get("probtables", Json::Value()).get("ins.len", Json::Value()).get("max.len", DEFAULT_MAX_INS_LENGTH).asUInt64() + 1));
            container->addDataVector(vector<prob_t>(_config.get("probtables", Json::Value()).get("ins.len", Json::Value()).get("max.len", DEFAULT_MAX_INS_LENGTH).asUInt64() + 1));
            containers[VDJ_VAR_DIV_INS_LEN] = container;

            // VD nuc + DJ nuc
            container = new TDVectorList(true, _config.get("probtables", Json::Value()).get("ins.nucl", Json::Value()).get("laplace", .0).asDouble());
            for (auto i = 0; i < 8; ++i) {
                container->addColumnName("VD/DJ nucs");
                container->addDataVector(vector<prob_t>(4));
            }
            containers[VDJ_VAR_DIV_INS_NUC] = container;
        }


        bool makeModelParameterVector(vector<AbstractTDContainer*> &containers) {
            // Made ModelParameterVector from input tables if all is ok.
            vector<prob_t> event_probs;  // param vec
            vector<event_ind_t> event_lengths;  // lens vec
            vector<event_ind_t> event_classes;  // event classes
            vector<seq_len_t> event_col_num;  // event family col numbers
            vector<prob_t> laplace;
            vector<seq_len_t> min_D_len_vec;

            bool is_ok = false;

            if (_recomb == VJ_RECOMB) {
                is_ok = this->makeVJModelParameterVector(containers, event_probs, event_lengths, event_classes, event_col_num, laplace);
            }
            else if (_recomb == VDJ_RECOMB) {
                is_ok = this->makeVDJModelParameterVector(containers, event_probs, event_lengths, event_classes, event_col_num, laplace, min_D_len_vec);
            }

            // Free all memory for containers.
            for (size_t i = 0; i < containers.size(); ++i) {
                if (containers[i]) { delete containers[i]; }
            }

            return is_ok;
        }


        bool makeVJModelParameterVector(vector<AbstractTDContainer*> &containers,
                                        vector<prob_t> &event_probs,
                                        vector<event_ind_t> &event_lengths,
                                        vector<event_ind_t> &event_classes,
                                        vector<seq_len_t> &event_col_num,
                                        vector<prob_t> &laplace) {
            bool is_ok = false;

            if (containers[VJ_VAR_JOI_GEN]
                && containers[VJ_VAR_DEL]
                && containers[VJ_JOI_DEL]
                && containers[VJ_VAR_JOI_INS_LEN]
                && containers[VJ_VAR_JOI_INS_NUC]) {

                this->addGenes(containers[VJ_VAR_JOI_GEN],
                               _genes->V(),
                               _genes->J(),
                               event_probs,
                               event_lengths,
                               event_classes,
                               event_col_num,
                               laplace,
                               0);

                this->addDels(containers[VJ_VAR_DEL],
                              _genes->V(),
                              event_probs,
                              event_lengths,
                              event_classes,
                              event_col_num,
                              laplace,
                              1);

                this->addDels(containers[VJ_JOI_DEL],
                              _genes->J(),
                              event_probs,
                              event_lengths,
                              event_classes,
                              event_col_num,
                              laplace,
                              containers[VJ_VAR_DEL]->n_columns());

                this->addIns(containers[VJ_VAR_JOI_INS_LEN],
                             event_probs,
                             event_lengths,
                             event_classes,
                             event_col_num,
                             laplace,
                             containers[VJ_JOI_DEL]->n_columns(),
                             _config.get("probtables", Json::Value()).get("ins.len", Json::Value()).get("max.len", DEFAULT_MAX_INS_LENGTH).asUInt64() + 1);

                this->addIns(containers[VJ_VAR_JOI_INS_NUC],
                             event_probs,
                             event_lengths,
                             event_classes,
                             event_col_num,
                             laplace,
                             1);

                _param_vec = new ModelParameterVector(VJ_RECOMB, event_probs, event_lengths, event_classes, event_col_num, laplace);
                is_ok = true;
            }

            return is_ok;
        }


        bool makeVDJModelParameterVector(vector<AbstractTDContainer*> &containers,
                                         vector<prob_t> &event_probs,
                                         vector<event_ind_t> &event_lengths,
                                         vector<event_ind_t> &event_classes,
                                         vector<seq_len_t> &event_col_num,
                                         vector<prob_t> &laplace,
                                         vector<seq_len_t> &min_D_len_vec) {
            bool is_ok = false;

            if (containers[VDJ_VAR_GEN]
                && containers[VDJ_JOI_DIV_GEN]
                && containers[VDJ_VAR_DEL]
                && containers[VDJ_JOI_DEL]
                && containers[VDJ_DIV_DEL]
                && containers[VDJ_VAR_DIV_INS_LEN]
                && containers[VDJ_VAR_DIV_INS_NUC]) {

                this->addGenes(containers[VDJ_VAR_GEN],
                               _genes->V(),
                               event_probs,
                               event_lengths,
                               event_classes,
                               event_col_num,
                               laplace);

                this->addGenes(containers[VDJ_JOI_DIV_GEN],
                               _genes->J(),
                               _genes->D(),
                               event_probs,
                               event_lengths,
                               event_classes,
                               event_col_num,
                               laplace,
                               1);

                this->addDels(containers[VDJ_VAR_DEL],
                              _genes->V(),
                              event_probs,
                              event_lengths,
                              event_classes,
                              event_col_num,
                              laplace,
                              1);

                this->addDels(containers[VDJ_JOI_DEL],
                              _genes->J(),
                              event_probs,
                              event_lengths,
                              event_classes,
                              event_col_num,
                              laplace,
                              containers[VDJ_VAR_DEL]->n_columns());

                this->addDels2D(containers[VDJ_DIV_DEL],
                                _genes->D(),
                                event_probs,
                                event_lengths,
                                event_classes,
                                event_col_num,
                                laplace,
                                containers[VDJ_JOI_DEL]->n_columns());

                this->addIns(containers[VDJ_VAR_DIV_INS_LEN],
                             event_probs,
                             event_lengths,
                             event_classes,
                             event_col_num,
                             laplace,
                             containers[VDJ_DIV_DEL]->n_rows(),
                             _config.get("probtables", Json::Value()).get("ins.len", Json::Value()).get("max.len", DEFAULT_MAX_INS_LENGTH).asUInt64() + 1);

                this->addIns(containers[VDJ_VAR_DIV_INS_NUC],
                             event_probs,
                             event_lengths,
                             event_classes,
                             event_col_num,
                             laplace,
                             1);

                for (seg_index_t i = 1; i <= _genes->D().max(); ++i) { min_D_len_vec.push_back(_min_D_len); }
                _param_vec = new ModelParameterVector(VDJ_RECOMB, event_probs, event_lengths, event_classes, event_col_num, laplace, true, min_D_len_vec);
                is_ok = true;
            }

            return is_ok;
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