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
        ProbabilisticAssemblingModel(const string& folderpath) : _model_path(folderpath + "/") {
            _vj_recomb = false;

            _status = false;
            _genes = nullptr;
            _param_vec = nullptr;

            _builder = nullptr;
            _generator = nullptr;

            _status = this->parseModelConfig(folderpath + "/model.json")
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
//            return this->_builder->buildAndCompute(???);
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
        MAAGRepertoire buildGraphs(const ClonesetView& repertoire,
                                   bool full_build = true,
                                   bool aminoacid = false) const {
            _builder->build(repertoire);
        }


        /**
         * \brief Generate artificial repertoire of sequences using this model.
         *
         * \param count Number of generated sequences.
         * \param merge If true than return only unique clonotypes with counts of occurrences for each clonotype.
         *
         * \return Artificial repertoire.
         */
        Cloneset generateSequences(size_t count = 1, bool merge = false) const {
            // generate sequences
//            return _generator->generate(count);
        }


        /**
         * \brief Access to the gene segments table.
         *
         * \return Const reference to the gene segments table.
         */
        const VDJRecombinationGenes& gene_segments() const { return *_genes; }


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

    protected:

        bool _status, _verbose;

        Json::Value _config;
        bool _vj_recomb;
        string _model_path;
        VDJRecombinationGenes *_genes;
        ModelParameterVector *_param_vec;

        MAAGBuilder *_builder;
        ClonotypeAssembler *_generator;


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
                        endl;
                return true;
            }
            cerr << "Probabilistic assembling model error:" << endl << "\t config .json file not found at [" << jsonpath << "]" << endl;
            return false;
        }


        /**
         * \brief Parse gene segment JSON files and tables.
         */
        bool parseGeneSegments() {
            if (_config.get("segments", Json::Value("")).get("variable", Json::Value("")).size() == 0) {
                cerr << "Probabilistic assembling model error:" << endl << "\t V(ariable) gene segments file not found" << endl;
                return false;
            }
            string v_path = _model_path + _config.get("segments", Json::Value("")).get("variable", Json::Value("")).get("file", "").asString();


            if (_config.get("segments", Json::Value("")).get("joining", Json::Value("")).size() == 0) {
                cerr << "Probabilistic assembling model error:" << endl << "\t J(oining) gene segments file not found" << endl;
                return false;
            }
            string j_path = _model_path + _config.get("segments", Json::Value("")).get("joining", Json::Value("")).get("file", "").asString();


            bool vok, jok, dok = true;
            if (_config.get("segments", Json::Value("")).get("diversity", Json::Value("")).size() == 0) {
                _genes = new VDJRecombinationGenes("VJ.V", v_path, "VJ.J", j_path, &vok, &jok);
                /* add P nucleotides */
            } else {
                string d_path = _model_path + _config.get("segments", Json::Value("")).get("diversity", Json::Value("")).get("file", "").asString();
                _genes = new VDJRecombinationGenes("VDJ.V", v_path, "VDJ.J", j_path, "VDJ.D", d_path, &vok, &jok, &dok);
                /* add P nucleotides */
            }

            return vok && jok && dok;
        }


        /**
         * \brief Parse files with event probabilities matrices and make ModelParameterVector.
         */
        bool parseEventProbabilities() {
            Json::Value pt = _config.get("probtables", "no-prob");
            if (pt.size()) {
                AbstractTDContainer *container;
                string element = "";
                for (Json::ArrayIndex i = 0; i < pt.size(); ++i) {
                    if (container) { delete container; }

                    element = pt.getMemberNames()[i];
                    container = read_textdata(pt[element]["file"].asString(),
                                  pt[element]["type"].asString(),
                                  pt[element]["skip.first.column"].asBool(),
                                  pt[element]["laplace"].asDouble());

                    if (container) {
                        if (_vj_recomb) {
                            if (element == "v.j") {

                            } else if (element == "v.del") {

                            } else if (element == "j.del") {

                            } else if (element == "ins.len") {

                            } else if (element == "ins.nucl") {

                            } else {
                                cerr << "Unrecognised element in \'probtables\'" << ":\n\t" << element << endl;
                            }
                        } else {
                            if (element == "v") {

                            } else if (element == "j.d") {

                            } else if (element == "v.del") {

                            } else if (element == "j.del") {

                            } else if (element == "d.del") {

                            } else if (element == "ins.len") {

                            } else if (element == "ins.nucl") {

                            } else {
                                cerr << "Unrecognised element in \'probtables\'" << ":\n\t" << element << endl;
                            }
                        }
                    }
                }

                if (container) { delete container; }

            } else {
                cerr << "No information about probability events in the model .json file found." << endl;
                return false;
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
            /*_generator = new ClonotypeAssembler(*_param_vec, *_genes);*/
        }

    };

}


#endif