//
// Created by Vadim N. on 01/03/2016.
//

#ifndef YMIR_MODEL_PARSER_H
#define YMIR_MODEL_PARSER_H


#include "textdata.h"
#include "genesegment.h"


namespace ymir {

    class ModelParser {

    public:

        ModelParser(const std::string &model_path, json_value config, ModelBehaviour behav)
            : _model_path(model_path),
              _config(config),
              _behaviour(behav)
        {
        }


        virtual ~ModelParser()
        {
        }


        bool parse() {
            return this->parseGeneSegments()
                   && this->makeEventProbabilitiesVector();
        }


        virtual bool parseGeneSegments() = 0;


        bool makeEventProbabilitiesVector() {
            if (_behaviour == EMPTY) {
                return this->createEventProbabilitiesFromScratch();
            } else {
                return this->parseEventProbabilitiesFromFiles();
            }
        }


        void swap_genes(unique_ptr<VDJRecombinationGenes> &ptr) {
            _genes.swap(ptr);
        }

        void swap_parameters(unique_ptr<ModelParameterVector> &ptr) {
            _param_vec.swap(ptr);
        }

    protected:

        json_value _config;
        std::string _model_path;
        ModelBehaviour _behaviour;
        unique_ptr<VDJRecombinationGenes> _genes;
        unique_ptr<ModelParameterVector> _param_vec;
        seq_len_t _min_D_len;


        bool createEventProbabilitiesFromScratch() {
            vector<AbstractTDContainer*> containers;
            containers.resize(10, nullptr);

            this->createContainers(containers);

            bool is_ok = this->makeModelParameterVector(containers);

            _param_vec->fill(1);
            _param_vec->normaliseEventFamilies();

            if (!is_ok) { throw(runtime_error("WRONG EMPTY VECTOR CREATING SUBROUTINE!!!")); }

            return is_ok;
        }


        virtual void createContainers(vector<AbstractTDContainer*> &containers) = 0;


        bool parseEventProbabilitiesFromFiles() {
            if (_config.find("probtables") != _config.end()) {
                json_value pt = _config["probtables"];

                AbstractTDContainer *container = nullptr;
                string element = "", err_message = "";
                vector<AbstractTDContainer*> containers;
                containers.resize(10, nullptr);

                // Parse tables with event probabilities.
                for (size_t i = 0; i < pt.size(); ++i) {
                    for (json::iterator it = pt.begin(); it != pt.end(); ++it) {
                        element = it.key();
                        container = read_textdata(_model_path + pt[element]["file"].get<std::string>(),
                                                  pt[element]["type"].get<std::string>(),
                                                  pt[element].value("skip.first.column", true),
                                                  pt[element].value("laplace", .0),
                                                  err_message);

                        this->parseDataContainer(element, container, containers);
                    }
                }

                return this->makeModelParameterVector(containers);
            }

            std::cout << "[ERROR] No information about probability events in the model .json file found." << std::endl;
            return false;
        }


        virtual void parseDataContainer(const string &element, AbstractTDContainer *container, vector<AbstractTDContainer*> &containers) = 0;


        bool makeModelParameterVector(vector<AbstractTDContainer*> &containers) {
            // Made ModelParameterVector from input tables if all is ok.
            vector<prob_t> event_probs;  // param vec
            vector<event_ind_t> event_lengths;  // lens vec
            vector<event_ind_t> event_classes;  // event classes
            vector<seq_len_t> event_col_num;  // event family col numbers
            vector<prob_t> laplace;
            vector<seq_len_t> min_D_len_vec;

            bool is_ok = this->makeModelParameterVector(containers, event_probs, event_lengths, event_classes, event_col_num, laplace, min_D_len_vec);

            // Free all memory for containers.
            for (size_t i = 0; i < containers.size(); ++i) {
                if (containers[i]) { delete containers[i]; }
            }

            return is_ok;
        }


        virtual bool makeModelParameterVector(vector<AbstractTDContainer*> &containers,
                                              vector<prob_t> &event_probs,
                                              vector<event_ind_t> &event_lengths,
                                              vector<event_ind_t> &event_classes,
                                              vector<seq_len_t> &event_col_num,
                                              vector<prob_t> &laplace,
                                              vector<seq_len_t> &min_D_len_vec) = 0;


        /**
         *
         */
        bool findGenes(const vector<string> &names, const GeneSegmentAlphabet &gsa, string &err_message) const {
            unordered_set<string> nameset;

            for (size_t i = 0; i < names.size(); ++i) {
                if (gsa[names[i]].index == 0) {
                    err_message = "ERROR: can't find '" + names[i] + "' in gene segments.";
                    return false;
                }
                nameset.insert(names[i]);
            }

            for (size_t i = 1; i <= gsa.max(); ++i) {
                if (nameset.count(gsa[i].allele) == 0) {
                    err_message = "ERROR: can't find '" + gsa[i].allele + "' in this file.";
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

                // TODO: what the hell is this!? Why do I need this?!
//                if (prob_data.size() > gsa[name_order[j]].sequence.size() + 1) {
//                    vector<prob_t> new_prob_data((gsa[name_order[j]].sequence.size() + 1) * (gsa[name_order[j]].sequence.size() + 1));
//                    for (auto row_i = 0; row_i < gsa[name_order[j]].sequence.size() + 1; ++row_i) {
//                        for (auto col_i = 0; col_i < gsa[name_order[j]].sequence.size() + 1; ++col_i) {
//                            new_prob_data[row_i * (gsa[name_order[j]].sequence.size() + 1) + col_i] = prob_data[row_i * (gsa[name_order[j]].sequence.size() + 1) + col_i];
//                        }
//                    }
//                    prob_data = new_prob_data;
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

    };


    class VJModelParser : public ModelParser {
    public:

        VJModelParser(const std::string &model_path, json_value config, ModelBehaviour behav)
            : ModelParser(model_path, config, behav)
        {
        }


        bool parseGeneSegments() {
            cout << "\tV gene seg.:     ";
            if (_config.value("segments", json_value("")).value("variable", json_value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                cout << "OK" << endl;
            }
            string v_path = _model_path + _config.value("segments", json_value("")).value("variable", json_value("")).value("file", "");

            cout << "\tJ gene seg.:     ";
            if (_config.value("segments", json_value("")).value("joining", json_value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                cout << "OK" << endl;
            }
            string j_path = _model_path + _config.value("segments", json_value("")).value("joining", json_value("")).value("file", "");

            bool vok, jok;

            _genes.reset(new VDJRecombinationGenes("VJ.V", v_path, "VJ.J", j_path, &vok, &jok));

            if (vok && jok) {
                _genes->appendPalindromicNucleotides(VARIABLE,
                                                     _config.value("segments", json_value("")).value("variable", json_value("")).value("P.nuc.3'", 0),
                                                     _config.value("segments", json_value("")).value("variable", json_value("")).value("P.nuc.5'", 0));
                _genes->appendPalindromicNucleotides(JOINING,
                                                     _config.value("segments", json_value("")).value("joining", json_value("")).value("P.nuc.3'", 0),
                                                     _config.value("segments", json_value("")).value("joining", json_value("")).value("P.nuc.5'", 0));
            }

            return vok && jok;
        }

    protected:

        virtual void createContainers(vector<AbstractTDContainer*> &containers) {
            AbstractTDContainer* container;

            // V-J
            container = new TDMatrix(true, _config.value("probtables", json_value("")).value("v.j", json_value("")).value("laplace", .0));
            for (auto i = 1; i <= _genes->V().max(); ++i){
                container->addRowName(_genes->V()[i].allele);
            }
            for (auto i = 1; i <= _genes->J().max(); ++i){
                container->addColumnName(_genes->J()[i].allele);
            }
            container->addDataVector(vector<prob_t>(_genes->V().max() * _genes->J().max()));
            containers[VJ_VAR_JOI_GEN] = container;

            // V del
            container = new TDVectorList(true, _config.value("probtables", json_value("")).value("v.del", json_value("")).value("laplace", .0));
            for (auto i = 1; i <= _genes->V().max(); ++i) {
                container->addColumnName(_genes->V()[i].allele);
                container->addDataVector(vector<prob_t>(_genes->V()[i].sequence.size() + 1));
            }
            containers[VJ_VAR_DEL] = container;

            // J del
            container = new TDVectorList(true, _config.value("probtables", json_value("")).value("j.del", json_value("")).value("laplace", .0));
            for (auto i = 1; i <= _genes->J().max(); ++i){
                container->addColumnName(_genes->J()[i].allele);
                container->addDataVector(vector<prob_t>(_genes->J()[i].sequence.size() + 1));
            }
            containers[VJ_JOI_DEL] = container;

            // VJ ins
            container = new TDVectorList(true, _config.value("probtables", json_value("")).value("ins.len", json_value("")).value("laplace", .0));
            container->addColumnName("VJ ins len");
            container->addDataVector(vector<prob_t>(_config.value("probtables", json_value("")).value("ins.len", json_value("")).value("max.len", DEFAULT_MAX_INS_LENGTH) + 1));
            containers[VJ_VAR_JOI_INS_LEN] = container;

            // VJ nuc
            container = new TDVectorList(true, _config.value("probtables", json_value("")).value("ins.nucl", json_value("")).value("laplace", .0));
            container->addColumnName("VJ nucs");
            container->addDataVector(vector<prob_t>(4));
            containers[VJ_VAR_JOI_INS_NUC] = container;
        }


        virtual void parseDataContainer(const string &element, AbstractTDContainer *container, vector<AbstractTDContainer*> &containers) {
            std::string err_message = "NOT FOUND";
            if (element == "v.j") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->J(), err_message)
                    && this->findGenes(container->row_names(), _genes->V(), err_message))
                {
                    containers[VJ_VAR_JOI_GEN] = container;
                    err_message = "OK";
                }

                cout << "\tV-J gene pairs:  " << err_message << endl;
            }
            else if (element == "v.del") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->V(), err_message))
                {
                    containers[VJ_VAR_DEL] = container;
                    err_message = "OK";
                }

                cout << "\tV delet. num.:   " << err_message << endl;
            }
            else if (element == "j.del") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->J(), err_message))
                {
                    containers[VJ_JOI_DEL] = container;
                    err_message = "OK";
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
                        err_message = "OK";
                    }
                }

                cout << "\tVJ ins. len.:    " << err_message << endl;
            }
            else if (element == "ins.nucl") {
                if (container && container->file_exists()) {
                    if (container->n_rows() != 4 || container->n_columns() != 1) {
                        stringstream ss;
                        ss << "ERROR: wrong number of columns and rows (expected: 4 X 1, got: " << (int) container->n_rows() << " X " << (int) container->n_columns() << ")";
                        err_message = ss.str();
                    } else {
                        containers[VJ_VAR_JOI_INS_NUC] = container;
                        err_message = "OK";
                    }
                }

                cout << "\tVJ ins. nuc.:    " << err_message << endl;
            }
            else { cerr << "Unrecognised element in \'probtables\'" << ":\n\t" << element << endl; }
        }


        virtual bool makeModelParameterVector(vector<AbstractTDContainer*> &containers,
                                              vector<prob_t> &event_probs,
                                              vector<event_ind_t> &event_lengths,
                                              vector<event_ind_t> &event_classes,
                                              vector<seq_len_t> &event_col_num,
                                              vector<prob_t> &laplace,
                                              vector<seq_len_t> &min_D_len_vec)
        {
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
                             _config.value("probtables", json_value("")).value("ins.len", json_value("")).value("max.len", DEFAULT_MAX_INS_LENGTH) + 1);

                this->addIns(containers[VJ_VAR_JOI_INS_NUC],
                             event_probs,
                             event_lengths,
                             event_classes,
                             event_col_num,
                             laplace,
                             1);

                _param_vec.reset(new ModelParameterVector(VJ_RECOMB, event_probs, event_lengths, event_classes, event_col_num, laplace, _config.value("errors", 0)));
                is_ok = true;
            }

            return is_ok;
        }

    };


    class VDJModelParser : public ModelParser {

    public:

        VDJModelParser(const std::string &model_path, json_value config, ModelBehaviour behav)
                : ModelParser(model_path, config, behav)
        {
        }


        bool parseGeneSegments() {
            cout << "\tV gene seg.:     ";
            if (_config.value("segments", json_value("")).value("variable", json_value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                cout << "OK" << endl;
            }
            string v_path = _model_path + _config.value("segments", json_value("")).value("variable", json_value("")).value("file", "");

            cout << "\tJ gene seg.:     ";
            if (_config.value("segments", json_value("")).value("joining", json_value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                cout << "OK" << endl;
            }
            string j_path = _model_path + _config.value("segments", json_value("")).value("joining", json_value("")).value("file", "");

            bool vok, jok, dok = true;

            cout << "\tD gene seg.:     ";
            if (_config.value("segments", json_value("")).value("diversity", json_value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                string d_path = _model_path + _config.value("segments", json_value("")).value("diversity", json_value("")).value("file", "");
                _genes.reset(new VDJRecombinationGenes("VDJ.V", v_path, "VDJ.J", j_path, "VDJ.D", d_path, &vok, &jok, &dok));
                _min_D_len = _config.value("segments", json_value("")).value("diversity", json_value("")).value("min.len", DEFAULT_DIV_GENE_MIN_LEN);
                cout << "OK" << endl;
            }

            if (vok && jok && dok) {
                _genes->appendPalindromicNucleotides(VARIABLE,
                                                     _config.value("segments", json_value("")).value("variable", json_value("")).value("P.nuc.3'", 0),
                                                     _config.value("segments", json_value("")).value("variable", json_value("")).value("P.nuc.5'", 0));
                _genes->appendPalindromicNucleotides(JOINING,
                                                     _config.value("segments", json_value("")).value("joining", json_value("")).value("P.nuc.3'", 0),
                                                     _config.value("segments", json_value("")).value("joining", json_value("")).value("P.nuc.5'", 0));
                if (_genes->is_vdj()) {
                    _genes->appendPalindromicNucleotides(DIVERSITY,
                                                         _config.value("segments", json_value("")).value("diversity", json_value("")).value("P.nuc.3'", 0),
                                                         _config.value("segments", json_value("")).value("diversity", json_value("")).value("P.nuc.5'", 0));
                }
            }

            return vok && jok && dok;
        }

    protected:

        virtual void createContainers(vector<AbstractTDContainer*> &containers) {
            AbstractTDContainer* container;

            // V
            container = new TDVector(true, _config.value("probtables", json_value("")).value("v", json_value("")).value("laplace", .0));
            container->addDataVector(vector<prob_t>());
            for (auto i = 1; i <= _genes->V().max(); ++i) {
                container->addRowName(_genes->V()[i].allele);
                container->addDataValue(1);
            }
            containers[VDJ_VAR_GEN] = container;
            cout << "\tV genes prob.:   " << "CREATED" << endl;

            // J-D
            container = new TDMatrix(true, _config.value("probtables", json_value("")).value("j.d", json_value("")).value("laplace", .0));
            for (auto i = 1; i <= _genes->J().max(); ++i) {
                container->addRowName(_genes->J()[i].allele);
            }
            for (auto i = 1; i <= _genes->D().max(); ++i) {
                container->addColumnName(_genes->D()[i].allele);
            }
            container->addDataVector(vector<prob_t>(_genes->D().max() * _genes->J().max()));
            containers[VDJ_JOI_DIV_GEN] = container;
            cout << "\tJ-D gene pairs:  " << "CREATED" << endl;

            // V del
            container = new TDVectorList(true, _config.value("probtables", json_value("")).value("v.del", json_value("")).value("laplace", .0));
            for (auto i = 1; i <= _genes->V().max(); ++i) {
                container->addColumnName(_genes->V()[i].allele);
                container->addDataVector(vector<prob_t>(_genes->V()[i].sequence.size() + 1));
            }
            containers[VDJ_VAR_DEL] = container;
            cout << "\tV delet. num.:   " << "CREATED" << endl;;

            // J del
            container = new TDVectorList(true, _config.value("probtables", json_value("")).value("j.del", json_value("")).value("laplace", .0));
            for (auto i = 1; i <= _genes->J().max(); ++i) {
                container->addColumnName(_genes->J()[i].allele);
                container->addDataVector(vector<prob_t>(_genes->J()[i].sequence.size() + 1));
            }
            containers[VDJ_JOI_DEL] = container;
            cout << "\tJ delet. num.:   " << "CREATED" << endl;

            // D del
            container = new TDMatrixList(true, _config.value("probtables", json_value("")).value("d.del", json_value("")).value("laplace", .0));
            for (auto i = 1; i <= _genes->D().max(); ++i) {
                container->addColumnName(_genes->D()[i].allele);
                container->addDataVector(vector<prob_t>( (_genes->D()[i].sequence.size() + 1) * (_genes->D()[i].sequence.size() + 1) ));
                container->addRowName(_genes->D()[i].allele);
                container->addMetadata(_genes->D()[i].sequence.size() + 1);
            }
            containers[VDJ_DIV_DEL] = container;
            cout << "\tD delet. num.:   " << "CREATED" << endl;

            // VD ins + DJ ins
            container = new TDVectorList(true, _config.value("probtables", json_value("")).value("ins.len", json_value("")).value("laplace", .0));
            container->addColumnName("VD ins");
            container->addColumnName("DJ ins");
            container->addDataVector(vector<prob_t>(_config.value("probtables", json_value("")).value("ins.len", json_value("")).value("max.len", DEFAULT_MAX_INS_LENGTH) + 1));
            container->addDataVector(vector<prob_t>(_config.value("probtables", json_value("")).value("ins.len", json_value("")).value("max.len", DEFAULT_MAX_INS_LENGTH) + 1));
            containers[VDJ_VAR_DIV_INS_LEN] = container;
            cout << "\tVD/DJ ins. len.: " << "CREATED" << endl;

            // VD nuc + DJ nuc
            container = new TDVectorList(true, _config.value("probtables", json_value("")).value("ins.nucl", json_value("")).value("laplace", .0));
            for (auto i = 0; i < 8; ++i) {
                container->addColumnName("VD/DJ nucs");
                container->addDataVector(vector<prob_t>(4));
            }
            containers[VDJ_VAR_DIV_INS_NUC] = container;
            cout << "\tVD/DJ ins. nuc.: " << "CREATED" << endl;
        }


        virtual void parseDataContainer(const string &element, AbstractTDContainer *container, vector<AbstractTDContainer*> &containers) {
            std::string err_message = "NOT FOUND";
            if (element == "v") {
                if (container && container->file_exists()) {
                    containers[VDJ_VAR_GEN] = container;
                    err_message = "OK";
                }
                cout << "\tV genes prob.:   " << err_message << endl;
            }
            else if (element == "j.d") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->D(), err_message)
                    && this->findGenes(container->row_names(), _genes->J(), err_message))
                {
                    containers[VDJ_JOI_DIV_GEN] = container;
                    err_message = "OK";
                }

                cout << "\tJ-D gene pairs:  " << err_message << endl;
            }
            else if (element == "v.del") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->V(), err_message))
                {
                    containers[VDJ_VAR_DEL] = container;
                    err_message = "OK";
                }

                cout << "\tV delet. num.:   " << err_message << endl;;
            }
            else if (element == "j.del") {
                if (container
                    && container->file_exists()
                    && this->findGenes(container->column_names(), _genes->J(), err_message))
                {
                    containers[VDJ_JOI_DEL] = container;
                    err_message = "OK";
                }

                cout << "\tJ delet. num.:   " << err_message << endl;
            }
            else if (element == "d.del") {
                if (container && container->file_exists()) {
                    containers[VDJ_DIV_DEL] = container;
                    err_message = "OK";
                }

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
                        err_message = "OK";
                    }
                }

                cout << "\tVD/DJ ins. len.: " << err_message << endl;
            }
            else if (element == "ins.nucl") {
                if (container && container->file_exists()) {
                    if (container->n_rows() != 4 || container->n_columns() != 8) {
                        stringstream ss;
                        ss << "ERROR: wrong number of columns and rows (expected: 4 X 8, got: " << (int) container->n_rows() << " X " << (int) container->n_columns() << ")";
                        err_message = ss.str();
                    } else {
                        containers[VDJ_VAR_DIV_INS_NUC] = container;
                        err_message = "OK";
                    }
                }

                cout << "\tVD/DJ ins. nuc.: " << err_message << endl;
            }
            else { std::cout << "Unrecognised element in \'probtables\'" << ":\n\t" << element << std::endl; }
        }


        virtual bool makeModelParameterVector(vector<AbstractTDContainer*> &containers,
                                              vector<prob_t> &event_probs,
                                              vector<event_ind_t> &event_lengths,
                                              vector<event_ind_t> &event_classes,
                                              vector<seq_len_t> &event_col_num,
                                              vector<prob_t> &laplace,
                                              vector<seq_len_t> &min_D_len_vec)
        {
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
                             _config.value("probtables", json_value("")).value("ins.len", json_value("")).value("max.len", DEFAULT_MAX_INS_LENGTH) + 1);

                this->addIns(containers[VDJ_VAR_DIV_INS_NUC],
                             event_probs,
                             event_lengths,
                             event_classes,
                             event_col_num,
                             laplace,
                             1);

                for (seg_index_t i = 1; i <= _genes->D().max(); ++i) { min_D_len_vec.push_back(_min_D_len); }
                _param_vec.reset(new ModelParameterVector(VDJ_RECOMB, event_probs, event_lengths, event_classes, event_col_num, laplace, _config.value("errors", 0), true, min_D_len_vec));
                is_ok = true;
            }

            return is_ok;
        }

    };


//    class VD2JModelParser : public ModelParser {
//
//    };
}


#endif //YMIR_MODEL_PARSER_H
