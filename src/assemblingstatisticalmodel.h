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


#include "maagbuilder.h"
#include "repertoire.h"
#include "genesegment.h"
#include "clonotypeassembler.h"


namespace ymir {

    /**
    * \function read_matrix
    *
    * \brief
    *
    * \param filepath Path to the file with matrix, separated by spaces / tabs.
    *
    * \return New struct with matrix.
    */
    NamedVectorArray read_matrix(const string& filepath) {
        vector<string> colnames, rownames;
        vector<prob_t> values;

        ifstream ifs;
        ifs.open(filepath);

        if (ifs.is_open()) {
            stringstream line_stream;
            string line, word;
            int col = 0;
            bool read_header = true;
            while (!ifs.eof()) {
                getline(ifs, line);
                if (line[0] != '\n') {
                    line_stream.str(line);
                    if (read_header) {
                        while (!line_stream.eof()) {
                            getline(line_stream, word, '\t');
                            colnames.push_back(word);
                        }
                        read_header = false;
                    } else {
                        getline(line_stream, word, '\t');
                        rownames.push_back(word);
                        while (!line_stream.eof()) {
                            getline(line_stream, word, '\t');

                        }
                        ++col;
                    }
                    line_stream.clear();
                }
            }
        } else {
            cerr << "Matrix parsing error:" << endl << "\tinput file [" << filepath << "] not found" << endl;
            return NamedVectorArray();
        }



        // read first line with names

        // put first element from each row to column names

        // remove cell at the intersection of the first row and the first column

        ifs.close();
//        return NamedMatrix(colnames, rownames, Matrix<prob_t, Dynamic, Dynamic, ColMajor>(values.data()));
    }


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
    * \class AssemblingStatisticalModel
    */
    class AssemblingStatisticalModel {

    public:

        /**
        *
        */
        AssemblingStatisticalModel(const string& folderpath) {
            _status = false;
            _genes = nullptr;
            _param_vec = nullptr;

            _builder = nullptr;
            _generator = nullptr;

//            _status = this->parseModelConfig(folderpath + "/model.json");
//            if (_status) {
//                _status = this->parseGeneSegments();
//                if (_status) {
//                    _status = this->parseEventProbabilities();
//                }
//            }

            _status = this->parseModelConfig(folderpath + "/model.json")
                      && this->parseGeneSegments()
                      && this->parseEventProbabilities();

            if (_status) {
                this->make_builder();
                this->make_assembler();
            }
        }


        /**
        *
        */
        virtual ~AssemblingStatisticalModel() {
            if (_genes) { delete _genes; }
            if (_param_vec) { delete _param_vec; }
            if (_builder) { delete _builder; }
            if (_generator) { delete _generator; }
        }


        /**
        * \brief Compute full assembling probabilities of clones from the given repertoire.
        *
        * \param repertoire Repertoire with alignments and event indices.
        * \param aminoacid What probabilities will be computed - of nucleotide or amino acid sequences. If
        * this parameter is false and some clone will have amino acid sequence, than error message will be generated to the output
        * and probability of this clone will set to zero.
        *
        * \return Vector of full assembling probabilities.
        */
        vector<numeric> computeFullProbabilities(const ClonesetView& repertoire, bool aminoacid = false) const {

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
        MAAGRepertoire buildGraphs(const ClonesetView& repertoire, bool full_build = true, bool aminoacid = false) const {
            _builder->build(repertoire);
        }


        /**
         * \brief Generate artificial repertoire of sequences using this model.
         *
         * \param count Number of generated sequences.
         * \param merge If true than return only unique clonotypes.
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
        VDJRecombinationGenes *_genes;
        ModelParameterVector *_param_vec;

        MAAGBuilder *_builder;
        ClonotypeAssembler *_generator;


        /**
         * \brief Private default constructor.
         */
        AssemblingStatisticalModel() {}


        /**
         * \brief Parse JSON file with model parameters.
         */
        virtual bool parseModelConfig(const string& jsonpath) {
            ifstream ifs;
            ifs.open(jsonpath);
            if (ifs.is_open()) {
                ifs >> _config;
                return true;
            }
            cerr << "Assembling statistical model error:" << endl << "\t config .json file not found" << endl;
            return false;
        }


        /**
         * \brief Parse gene segment JSON files and tables.
         */
        virtual bool parseGeneSegments() {
            string v_path = _config.get("segments", "genes").get("variable", "").asString(),
                   j_path = _config.get("segments", "genes").get("joining", "").asString(),
                   d_path = _config.get("segments", "genes").get("diversity", "").asString();

            if (v_path == "") {
                cerr << "Assembling statistical model error:" << endl << "\t V(ariable) gene segments file not found" << endl;
                return false;
            }

            if (j_path == "") {
                cerr << "Assembling statistical model error:" << endl << "\t J(oining) gene segments file not found" << endl;
                return false;
            }

            bool vok, jok, dok = true;
            if (d_path == "") {
                _genes = new VDJRecombinationGenes("VJ.V", v_path, "VJ.J", j_path, &vok, &jok);
            } else {
                _genes = new VDJRecombinationGenes("VDJ.V", v_path, "VDJ.J", j_path, "VDJ.D", d_path, &vok, &jok, &dok);
            }

            return vok && jok && dok;
        }


        /**
         * \brief Parse files with event probabilities matrices and make ModelParameterVector.
         */
        virtual bool parseEventProbabilities() {
            // for each vector check size of gene segments and remove trailing zeros
            // or get a message when some column is shorter than length

            cerr << "Some problems here" << endl;
        }


        // resize deletions vectors to gene segment length
        void removeTrailingZeros() const {}


        // grow deletions vectors (fill with zeros added rows) so
        // all vectors will be of same length
        void addTrailingZeros() const {}


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