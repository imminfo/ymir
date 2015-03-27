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
    * \function read_event_matrix
    *
    * \brief
    *
    * \param filepath Path to the file with matrix, separated by spaces / tabs.
    *
    * \return New struct with matrix with event probabilities and metadata.
    */
//    event_matrix_t read_event_matrix(const string& filepath) {
//
//    }





    // If input is layout:
    // input/parser -> layout -> model/builder -> graph
    // segment sequences required only for those segments w/o layout

    // If input is sequences:
    // input/parser -> sequences -> layout -> model/builder -> graph
    // inner aligner will do the work

    // <!!!> input for ymir is default (layouts or sequences), but various parsers in Python, maybe..?
    // therefore no pYmir package.
    // for rYmir - parsers in R, sadly, but input are data frames, which will be
    // transformed to Ymir's format by Rcpp.

    // therefore no need for template model
    // only abstract class and two childs for VJ- and VDJ-recombinations.


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


    //
    //  MODEL BUILDER!!!!
    //

    /**
    * \class AssemblingStatisticalModel
    */
    class AssemblingStatisticalModel : protected MAAG {

    public:

        /**
        *
        */
        AssemblingStatisticalModel(const string& folderpath) : MAAG(false) {
            _status = false;
            _genes = nullptr;
            _model_parameters = nullptr;

            _generator = nullptr;

            _status = this->parseModelConfig(folderpath + "/model.json");
            if (_status) {
                _status = this->parseGeneSegments();
                if (_status) {
                    _status = this->parseEventProbabilities();
                }
            }
        }


        /**
        *
        */
        virtual ~AssemblingStatisticalModel() {
            if (_genes) {
                delete _genes;
            }
        }


        /**
        * \brief Compute full assembling probabilities of clones from the given repertoire.
        *
        * \param repertoire Repertoire with alignments and event indices.
        * \param aminoacid What probabilities will be computed - for nucleotide or amino acid sequences. If
        * this parameter is false and some clone will have amino acid sequence, than error message will be generated to the output
        * and probability of this clone will set to zero.
        *
        * \return Vector of full assembling probabilities.
        */
        const vector<numeric>& computeFullProbabilities(const ClonesetView& repertoire, bool aminoacid = false) const {
            // COMPUTE ALL GRAPHS -> COMPUTE ALL PROBABILITIES
            vector<AssemblyGraph> graph_vec = this->buildGraphs(repertoire);
//            vector<AssemblyGraph> graph_vec(this->buildGraphs(repertoire));  // which is slower???
            vector<numeric> prob_vec;
            prob_vec.reserve(graph_vec.size());

            // OMP cycle
            for (size_t i = 0; i < graph_vec.size(); i++) {
                prob_vec[i] = graph_vec[i].fullProbability();
            }


            // FOR EACH: COMPUTE ONE GRAPH -> GET PROBABILITY
            // in cycle compute graph and get full probabilties? is it slower?
            // omp cycle
        }


        /**
        * \brief Build assembly graphs for the given repertoire with sequences alignments.
        *
        * \param repertoire Repertoire with alignments and event indices.
        *
        * \return Vector of assembly graphs.
        */
//        const vector<AssemblyGraph>& buildGraphs(const ClonalRepertoireView& repertoire, bool aminoacid = false) const {
//
//        }


        /**
        * \brief Generate artificial repertoire of sequences using this model.
        *
        * \return Artificial repertoire.
        */
        Cloneset generateSequences(size_t count = 1) const {
            if (!_generator) {
                // make generator
            }

            // generate sequences
//            return this->_generator->generate(count);
        }


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


        bool status() const { return _status; }

    protected:

        bool _status;

        Json::Value _config;
        VDJRecombinationGenes *_genes;
        ModelParameterVector *_model_parameters;

        ClonotypeAssembler *_generator;



        AssemblingStatisticalModel() : AssemblyGraph(false) {}


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
        * \brief Parser gene segment JSON files and tables.
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
        * \brief Parse files with event probabilities matrices.
        */
        virtual bool parseEventProbabilities() {
            // make indexer here
        }


        // resize deletions vectors to gene segment length
        void removeTrailingZeros() {}

        // grow deletions vectors (fill with zeros added rows) so
        // all vectors will be of same length
        void addTrailingZeros() {}

    };


    /*
    Builder:

    nucleotide sequence builder

    amino acid sequence builder

    build()
    Clone, Metadata, ModelParameterVector
    -> AssemblyGraph
    ClonalRepertoireView, Metadata, ModelParameterVector
    -> AssemblyGraphRepertoire

    replaceEventProbabilities()
    &AssemblyGraph, ModelParameterVector
    -> bool
    &AssemblyGraphRepertoire, ModelParameterVector
    -> bool
    */
    // AssemblyGraphBuilder
}


#endif