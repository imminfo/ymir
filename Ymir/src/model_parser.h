//
// Created by Vadim N. on 01/03/2016.
//

#ifndef YMIR_MODEL_PARSER_H
#define YMIR_MODEL_PARSER_H


namespace ymir {

    class ModelParser {

    public:

        ModelParser(const std::string &model_path, Json::Value config, ModelBehaviour behav)
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

        Json::Value _config;
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


        bool makeModelParameterVector(vector<AbstractTDContainer*> &containers) {
            // Made ModelParameterVector from input tables if all is ok.
            vector<prob_t> event_probs;  // param vec
            vector<event_ind_t> event_lengths;  // lens vec
            vector<event_ind_t> event_classes;  // event classes
            vector<seq_len_t> event_col_num;  // event family col numbers
            vector<prob_t> laplace;
            vector<seq_len_t> min_D_len_vec;

            bool is_ok = this->makeModelParameterVector(containers, event_probs, event_lengths, event_classes, event_col_num, laplace);

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
                                              vector<prob_t> &laplace) = 0;


        bool parseEventProbabilitiesFromFiles() = 0;



        makeModelParameterVector

    };


    class VJModelParser : public ModelParser {
    public:

        VJModelParser(const std::string &model_path, Json::Value config, ModelBehaviour behav)
            : ModelParser(model_path, config, behav)
        {
        }


        bool parseGeneSegments() {
            cout << "\tV gene seg.:     ";
            if (_config.get("segments", Json::Value("")).get("variable", Json::Value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                cout << "OK" << endl;
            }
            string v_path = _model_path + _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("file", "").asString();

            cout << "\tJ gene seg.:     ";
            if (_config.get("segments", Json::Value("")).get("joining", Json::Value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                cout << "OK" << endl;
            }
            string j_path = _model_path + _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("file", "").asString();

            bool vok, jok;

            _genes.reset(new VDJRecombinationGenes("VJ.V", v_path, "VJ.J", j_path, &vok, &jok));

            if (vok && jok) {
                _genes->appendPalindromicNucleotides(VARIABLE,
                                                     _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("P.nuc.3'", 0).asUInt(),
                                                     _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("P.nuc.5'", 0).asUInt());
                _genes->appendPalindromicNucleotides(JOINING,
                                                     _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("P.nuc.3'", 0).asUInt(),
                                                     _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("P.nuc.5'", 0).asUInt());
            }

            return vok && jok;
        }

    protected:



    };


    class VDJModelParser : public ModelParser {

    public:

        VDJModelParser(const std::string &model_path, Json::Value config, ModelBehaviour behav)
                : ModelParser(model_path, config, behav)
        {
        }


        bool parseGeneSegments() {
            cout << "\tV gene seg.:     ";
            if (_config.get("segments", Json::Value("")).get("variable", Json::Value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                cout << "OK" << endl;
            }
            string v_path = _model_path + _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("file", "").asString();

            cout << "\tJ gene seg.:     ";
            if (_config.get("segments", Json::Value("")).get("joining", Json::Value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                cout << "OK" << endl;
            }
            string j_path = _model_path + _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("file", "").asString();

            bool vok, jok, dok = true;

            cout << "\tD gene seg.:     ";
            if (_config.get("segments", Json::Value("")).get("diversity", Json::Value("")).size() == 0) {
                cout << "ERROR: no gene segments file in the model's .json." << endl;

                return false;
            } else {
                string d_path = _model_path + _config.get("segments", Json::Value("")).get("diversity", Json::Value("")).get("file", "").asString();
                _genes.reset(new VDJRecombinationGenes("VDJ.V", v_path, "VDJ.J", j_path, "VDJ.D", d_path, &vok, &jok, &dok));
                _min_D_len = _config.get("segments", Json::Value("")).get("diversity", Json::Value("")).get("min.len", DEFAULT_DIV_GENE_MIN_LEN).asUInt();
                cout << "OK" << endl;
            }

            if (vok && jok && dok) {
                _genes->appendPalindromicNucleotides(VARIABLE,
                                                     _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("P.nuc.3'", 0).asUInt(),
                                                     _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("P.nuc.5'", 0).asUInt());
                _genes->appendPalindromicNucleotides(JOINING,
                                                     _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("P.nuc.3'", 0).asUInt(),
                                                     _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("P.nuc.5'", 0).asUInt());
                if (_genes->is_vdj()) {
                    _genes->appendPalindromicNucleotides(DIVERSITY,
                                                         _config.get("segments", Json::Value("")).get("diversity", Json::Value("")).get("P.nuc.3'", 0).asUInt(),
                                                         _config.get("segments", Json::Value("")).get("diversity", Json::Value("")).get("P.nuc.5'", 0).asUInt());
                }
            }

            return vok && jok && dok;
        }

    protected:

    };


//    class VD2JModelParser : public ModelParser {
//
//    };
}


#endif //YMIR_MODEL_PARSER_H
