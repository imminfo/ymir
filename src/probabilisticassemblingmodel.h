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

    // ymir probs --repertoire twb.txt --model TRB --parser mitcr.json
    // ymir probs -r twb.txt -m TRB -p mitcr.json

    // ymir inference --repertoire twa.txt --model TRA --n_iter 100 --eps 1e-4 --algorithm online-em --parser mitcr.json
    // ymir inference -r twa.txt -m TRA -n 100 -e 1e-4 -a online-em -p mitcr.json

    // ymir generate --count 100000 --model TRB --out generated.txt
    // ymir generate -c 100000 -m TRB -o generated.txt


    // parsers:
    // MiTCR, MiXCR, MiGEC
    // V.segments / V.end
    // J.segments / J.start
    /*
    column names:
    --vseg --jseg
    --vstart --vend
    --seq
     */

    // etc.

    // case with a few segments, not only one:
    // choose best / get all / get parameters for all
    // i.e., inference w/ or w/o inference of segments usage

    /**
    * \class ProbabilisticAssemblingModel
    */
    class ProbabilisticAssemblingModel {

    public:

        /**
        *
        */
        ProbabilisticAssemblingModel(const string& folderpath) {
            _model_path = folderpath + "/";

            _vj_recomb = false;
            _status = false;

            _genes = nullptr;
            _param_vec = nullptr;

            _builder = nullptr;
            _generator = nullptr;

            _status = this->parseModelConfig(_model_path + "model.json")
                      && this->parseGeneSegments()
                      && this->parseEventProbabilities();

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
                                                MAAG_COMPUTE_PROB_ACTION action = MAX_PROBABILITY) const {
            return this->_builder->buildAndCompute(repertoire, aminoacid, action);
        }


        /**
         * \brief Build a set of MAAGs from the given cloneset.
         *
         * \param repertoire Cloneset.
         * \param full_builder If true than build full MAAG (with additional matrices for event indices).
         * \param aminoacid If true than build graphs of aminoacid CDR3 sequences, traslating the given repertoire first.
         *
         * \return Set of MAAGs.
         */
        MAAGRepertoire buildGraphs(const ClonesetView &repertoire,
                                   bool full_build = true,
                                   bool aminoacid = false) const {
            return _builder->build(repertoire, full_build);
        }


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
        const VDJRecombinationGenes& gene_segments() const { return *_genes; }


        /**
         * \brief Access to vector of probabilities.
         */
        const ModelParameterVector& event_probabilities() const { return *_param_vec; }


        /**
        * \brief Read the model from the given folder.
        *
        * Search for model.json file.
        */
        bool read(const string& folderpath) {
            // check for folder and return false if something wrong.
            // read json file with model
//            print name
//            print full name
//            print comment
        }


        /**
        * \brief Save the model to the given folder.
        *
        * \param folderpath Path to the model folder.
        *
        * \return True if all ok.
        */
        bool write(const string& folderpath) const {
            // get JSON object and write it to model.json file
            // get all prob tables and write them to files (aggregate those with similar names to matrix.list)
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

        Json::Value _config;
        bool _vj_recomb;
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
                _vj_recomb = _config.get("recombination", "VJ").asString() == "VJ";
                cout << "Probabilistic assembling model:\n\t" <<
                        _config.get("name", "Nameless model").asString() <<
                        "\n\t(" <<
                        _config.get("comment", "").asString() <<
                        ")\n\t" <<
                        _config.get("recombination", "VJ").asString() <<
                        "-recombination  |  " <<
                        (_config.get("hypermutations", false).asBool() ? "Hypermutations" : "No hypermutations") <<
                        endl << "\tFiles:"<< endl;
                return true;
            }
            cerr << "Probabilistic assembling model error:" << endl << "\t config .json file not found at [" << jsonpath << "]" << endl;
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

            if (_vj_recomb) {
                _genes = new VDJRecombinationGenes("VJ.V", v_path, "VJ.J", j_path, &vok, &jok);
                /* add P nucleotides */
            } else {
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
            }

            return vok && jok && dok;
        }


        /**
         * \brief Parse files with event probabilities matrices and make ModelParameterVector.
         */
        bool parseEventProbabilities() {
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
                                              pt[element]["laplace"].asDouble(),
                                              err_message);

                    if (_vj_recomb) {
                        bool check = true;
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
                                if (container->row_names().size() != 4 || container->column_names().size() != 4) {
                                    stringstream ss;
                                    ss << "ERROR: wrong number of columns and rows (expected: 4 X 4, got: " << (int) container->row_names().size() << "X" << (int) container->column_names().size() << ")";
                                    err_message = ss.str();
                                } else {
                                    containers[VJ_VAR_JOI_INS_NUC] = container;
                                }
                            }

                            cout << "\tVJ ins. nuc.:    " << err_message << endl;
                        }
                        else { cerr << "Unrecognised element in \'probtables\'" << ":\n\t" << element << endl; }
                    }
                    else {
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
                                if (container->row_names().size() != 4 || container->column_names().size() != 8) {
                                    stringstream ss;
                                    ss << "ERROR: wrong number of columns and rows (expected: 4 X 8, got: " << (int) container->row_names().size() << "X" << (int) container->column_names().size() << ")";
                                    err_message = ss.str();
                                } else {
                                    containers[VDJ_VAR_DIV_INS_NUC] = container;
                                }
                            }

                            cout << "\tVD/DJ ins. nuc.: " << err_message << endl;
                        }
                        else { cerr << "Unrecognised element in \'probtables\'" << ":\n\t" << element << endl; }
                    }
                }

                // Made ModelParameterVector from input tables if all is ok.
                vector<prob_t> event_probs;  // param vec
                vector<eventind_t> event_lengths;  // lens vec
                vector<eventind_t> event_classes;  // event classes
                vector<seq_len_t> event_col_num;  // event family col numbers
                vector<prob_t> laplace;
                vector<segindex_t> name_order;
                vector<seq_len_t> min_D_len_vec;

                bool is_ok = false;
                if (_vj_recomb) {
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
                                     containers[VJ_JOI_DEL]->n_columns());

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
                } else {
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
                                     containers[VDJ_DIV_DEL]->row_names().size());

                        this->addIns(containers[VDJ_VAR_DIV_INS_NUC],
                                     event_probs,
                                     event_lengths,
                                     event_classes,
                                     event_col_num,
                                     laplace,
                                     1);

                        for (segindex_t i = 0; i < _genes->D().size(); ++i) { min_D_len_vec.push_back(_min_D_len); }
                        _param_vec = new ModelParameterVector(VDJ_RECOMB, event_probs, event_lengths, event_classes, event_col_num, laplace, true, min_D_len_vec);
                        is_ok = true;
                    }
                }

                // Free all memory for containers.
                for (size_t i = 0; i < containers.size(); ++i) {
                    if (containers[i]) { delete containers[i]; }
                }

                return is_ok;

            } else {
                cerr << "No information about probability events in the model .json file found." << endl;
                return false;
            }
        }


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


        vector<segindex_t> arrangeNames(const vector<string> &names, const GeneSegmentAlphabet &gsa) const {
            vector<segindex_t> res;
            res.resize(names.size(), 0);

            for (size_t i = 0; i < names.size(); ++i) { res[i] = gsa[names[i]].index; }

            return res;
        }


        void addGenes(AbstractTDContainer *container,
                      const GeneSegmentAlphabet &gsa,
                      vector<prob_t> &event_probs,
                      vector<eventind_t> &event_lengths,
                      vector<eventind_t> &event_classes,
                      vector<seq_len_t> &event_col_num,
                      vector<prob_t> &laplace) const {
            vector<segindex_t> name_order = this->arrangeNames(container->row_names(), gsa);
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
                      vector<eventind_t> &event_lengths,
                      vector<eventind_t> &event_classes,
                      vector<seq_len_t> &event_col_num,
                      vector<prob_t> &laplace,
                      segindex_t prev_class_size) const {
            vector<segindex_t> name_order_row = this->arrangeNames(container->row_names(), gsa_row);
            vector<segindex_t> name_order_column = this->arrangeNames(container->column_names(), gsa_column);
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
                     vector<eventind_t> &event_lengths,
                     vector<eventind_t> &event_classes,
                     vector<seq_len_t> &event_col_num,
                     vector<prob_t> &laplace,
                     segindex_t prev_class_size) const {
            vector<segindex_t> name_order = this->arrangeNames(container->column_names(), gsa);
            vector<prob_t> prob_data;
            for (size_t i = 0; i < name_order.size(); ++i) {
                // find correct segment for i-th position
                size_t j = 0;
                for (; i+1 != name_order[j] ; ++j) {}
                prob_data = container->data(j);
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
                       vector<eventind_t> &event_lengths,
                       vector<eventind_t> &event_classes,
                       vector<seq_len_t> &event_col_num,
                       vector<prob_t> &laplace,
                       segindex_t prev_class_size) const {
            vector<segindex_t> name_order = this->arrangeNames(container->row_names(), gsa);
            vector<prob_t> prob_data;
            for (size_t i = 0; i < name_order.size(); ++i) {
                // find correct segment for i-th position
                size_t j = 0;
                for (; i+1 != name_order[j] ; ++j) {}
                prob_data = container->data(j);

//                if (prob_data.size() > gsa[name_order[j]].sequence.size() + 1) {
//                    prob_data.resize(gsa[name_order[j]].sequence.size() + 1);
//                }

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
                    vector<eventind_t> &event_lengths,
                    vector<eventind_t> &event_classes,
                    vector<seq_len_t> &event_col_num,
                    vector<prob_t> &laplace,
                    segindex_t prev_class_size) const {
            vector<prob_t> prob_data;
            for (size_t i = 0; i < container->n_columns(); ++i) {
                prob_data = container->data(i);
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

    };

}


#endif