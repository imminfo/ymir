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

#ifndef _PARSER_H
#define _PARSER_H


#include "iostream"
#include "fstream"
#include "sstream"

#include "aligner.h"
#include "repertoire.h"
#include "genesegment.h"


namespace ymir {

    #define DEFAULT_BUFFER_SIZE 100000


    /**
    * \class RepertoireParser
    *
    * \brief Parser for text files with repertoire data. By default it's MiTCR parser.
    * To make new parsers inherit from this class and rewrite virtual private method
    * "parseRepertoire".
    */
    class RepertoireParser {
    public:

        /**
        * \enum ALIGNMENT_COLUMN_ACTION
        *
        * \brief Specify an action to perform with an alignment column:
        * either overwrite found alignments (OVERWRITE), make alignments (MAKE_IF_NOT_FOUND) of clones if no such column found
        * or skip (SKIP) this alignments (e.g., D(iversity) gene segment for TCR alpha-chains).
        */
        enum ALIGNMENT_COLUMN_ACTION {
            OVERWRITE,
            MAKE_IF_NOT_FOUND,
            SKIP
        };


        struct AlignmentColumnOptions {
            ALIGNMENT_COLUMN_ACTION align_V, align_J, align_D;

            AlignmentColumnOptions() {}


            AlignmentColumnOptions& setV(ALIGNMENT_COLUMN_ACTION action) {
                this->align_V = action;
                return *this;
            }

            AlignmentColumnOptions& setJ(ALIGNMENT_COLUMN_ACTION action) {
                this->align_J = action;
                return *this;
            }

            AlignmentColumnOptions& setD(ALIGNMENT_COLUMN_ACTION action) {
                this->align_D = action;
                return *this;
            }
        };


        /**
        * \typedef ParserConfig
        *
        * \brief Parameters of this parser: parser's name, names of the columns, separators, sequences and other.
        */
        typedef Json::Value ParserConfig;


        RepertoireParser() {
//            _config_is_loaded = false;
        }


        /**
        * \brief Parse text file with sequence alignment data and return constructed Cloneset.
        *
        * Parse all lines in the file and return a repertoire. If no alignments found or any of align_* parameters
        * is true, create alignment using input aligner from function loadFile().
        *
        * \param rep Pointer to clonal repertoire object to which data will be uploaded.
        * \param gene_segments Alphabets of gene segments.
        * \param filepath Path to a file with sequences.
        * \param aligner Aligner object for alignment of clones.
        * \param nuc_sequences Boolean - if true than generate clones with nucleotide sequences from column
        * with nucleotide sequences data. Otherwise generate clones with amino acid sequences from corresponding
        * column.
        * \param align_V What action to do with column with V alignments.
        * \param align_J What action to do with column with J alignments.
        * \param align_D What action to do with column with D alignments.
        *
        * \return True if found has been successfully processed, false otherwise.
        */
        bool parse(const string& filepath,
                   Cloneset *rep,
                   const VDJRecombinationGenes& gene_segments,
                   AlignmentColumnOptions opts = AlignmentColumnOptions().setV(MAKE_IF_NOT_FOUND).setJ(MAKE_IF_NOT_FOUND).setD(MAKE_IF_NOT_FOUND),
                   const AbstractAligner& aligner = NaiveNucleotideAligner(),
                   bool nuc_sequences = true) {

//            if (!_config_is_loaded) { return false; }

            vector<Clonotype> clonevec;
            clonevec.reserve(DEFAULT_REPERTOIRE_RESERVE_SIZE);
            ifstream ifs;
            ifs.open(filepath);
            bool res = false;
            if (ifs.is_open()) {
                std::cout << "Parsing input file:\t" << filepath << endl;
                res = this->parseRepertoire(filepath, ifs, clonevec, gene_segments, aligner, nuc_sequences, opts);
                if (res) { rep->swap(clonevec); }
            } else {
                std::cout << "Repertoire parser error:" << "\tinput file [" << filepath << "] not found" << endl;
                res = false;
            }

            ifs.close();
            return res;
        }


    protected:

//        ParserConfig _config;
//        bool _config_is_loaded;


        virtual bool parseRepertoire(const string &filename,
                                     ifstream& ifs,
                                     vector<Clonotype>& vec,
                                     const VDJRecombinationGenes& gene_segments,
                                     const AbstractAligner& aligner,
                                     bool nuc_sequences,
                                     AlignmentColumnOptions opts) {

            unordered_map<string, int> columns_id;

            string col_nuc_seq = "Nucleotide sequence";
//            columns_id[col_nuc_seq] = -1;
            int col_nuc_seq_id = 0;

            string col_aa_seq = "Amino acid sequence";
//            columns_id[col_aa_seq] = -1;
            int col_aa_seq_id = 1;

            string col_v_seg = "Variable";
//            columns_id[col_v_seg] = -1;
            int col_v_seg_id = 2;

            string col_j_seg = "Joining";
//            columns_id[col_j_seg] = -1;
            int col_j_seg_id = 4;

            string col_d_seg = "Diversity";
//            columns_id[col_d_seg] = -1;
            int col_d_seg_id = 3;

            string col_v_end = "V end";
//            columns_id[col_v_end] = -1;
            int col_v_end_id = 5;

            string col_j_start = "J start";
//            columns_id[col_j_start] = -1;
            int col_j_start_id = 8;

            string col_d_als = "D.alignments";
//            columns_id[col_d_als] = -1;
            int col_d_als_id = 6;

            char sep = '\t';
            char al_sep = ',';
//            bool based1 = _config.get("1-based", false).asBool();

            bool do_align_V = false, do_align_J = false, do_align_D = false;

            stringstream line_stream, word_stream;
            string line, word, sequence;

            int index = 0;
            bool header = true;

            vector<seg_index_t> vseg, jseg, dseg;
            string temp_str;

            ClonotypeBuilder clone_builder;

            int glob_index = 1, bad_index = 0;

            while (!ifs.eof()) {
                getline(ifs, line);
                if (line.size() > 2) {
                    line_stream.str(line);
                    if (header) {
                        // check main column
                        if (nuc_sequences) {
                            if (columns_id[col_nuc_seq] == -1) {
                                cout << "Repertoire parser error:" << endl << "\tcolumn with nucleotide sequences not found" << endl;
                                return false;
                            }
                        } else {
                            if (columns_id[col_aa_seq] == -1) {
                                cout << "Repertoire parser error:" << endl << "\tcolumn with amino acid sequences not found" << endl;
                                return false;
                            }
                        }

                        // check for alignments
                        if ((columns_id[col_v_end] == -1 && opts.align_V == MAKE_IF_NOT_FOUND) || opts.align_V == OVERWRITE) {
                            do_align_V = true;
                        }
                        if ((columns_id[col_j_start] == -1 && opts.align_J == MAKE_IF_NOT_FOUND) || opts.align_J == OVERWRITE) {
                            do_align_J = true;
                        }
                        if ((gene_segments.is_vdj()) && ((columns_id[col_d_als] == -1 && opts.align_D == MAKE_IF_NOT_FOUND) || opts.align_D == OVERWRITE)) {
                            do_align_D = true;
                        }

                        header = false;

                    } else {

                        // parse body and build clones from each line

                        index = 0;

                        vseg.clear();
                        jseg.clear();
                        dseg.clear();

                        while (!line_stream.eof()) {

                            getline(line_stream, word, sep);

                            if (index == col_nuc_seq_id) {
                                if (nuc_sequences) {
                                    clone_builder.setSequence(word);
                                    clone_builder.setNucleotideSeq();
                                    sequence = word;
                                }
                            } else if (index == col_aa_seq_id) {
                                if (!nuc_sequences) {
                                    clone_builder.setSequence(word);
                                    clone_builder.setAminoAcidSeq();
                                    sequence = word;
                                }
                            } else if (index == col_v_seg_id) {
                                parseWordSegment(word, al_sep, vseg, gene_segments.V(), glob_index);
                            } else if (index == col_v_end_id) {
                                parseWordAlignment(word, al_sep, vseg, gene_segments.V(), clone_builder, 'V');
                            } else if (index == col_j_seg_id) {
                                parseWordSegment(word, al_sep, jseg, gene_segments.J(), glob_index);
                            } else if (index == col_j_start_id) {
                                parseWordAlignment(word, al_sep, jseg, gene_segments.J(), clone_builder, 'J');
                            } else if (index == col_d_seg_id) {
                                // NOT IMPLEMENTED YET
                                if (!do_align_D) {
                                    // error if VDJ recombination, no D alignment and D genes column is bad
                                }
                            } else if (index == col_d_als_id) {
                                if (!do_align_D) {
                                    // error if VDJ recombination and D alignment columns is bad
                                }
                            }

                            ++index;
                        }

                        bool align_ok = false;
                        if (do_align_V) {}

                        if (do_align_J) {}

                        if (do_align_D) {
                            for (seg_index_t seg_i = 1; seg_i <= gene_segments.D().max(); ++seg_i) {
                                AbstractAligner::LocalAlignmentIndices indices =
                                        aligner.alignLocal(gene_segments.D()[seg_i].sequence,
                                                           sequence,
                                                           DEFAULT_DIV_GENE_MIN_LEN);
                                if (indices.size()) {
                                    align_ok = true;
                                    for (size_t align_i = 0; align_i < indices.size(); ++align_i) {
                                        clone_builder.addDalignment(seg_i, indices[align_i]);
                                    }
                                }
                            }
                            if (!align_ok) {
                                cout << "Diversity gene could NOT be aligned with the given minimal gene length (min gene length " << (size_t) DEFAULT_DIV_GENE_MIN_LEN << ", line " << (size_t) glob_index << ")" << endl;
                                ++bad_index;
                            }
                        }

                        if (glob_index % 50000 == 0) {
                            cout << this->get_prefix(filename) + "parsed " << (size_t) glob_index << " lines" << endl;
                        }
                        ++glob_index;

                        //
                        // remove bad clonotypes here
                        //

//                        if (align_ok) {
                        vec.push_back(clone_builder.buildClonotype());
//                        }
                    }
                    line_stream.clear();
                }
            }

            std::cout << this->get_prefix(filename) + "parsed " <<
                    (size_t) glob_index <<
                    " lines (" <<
                    (size_t) bad_index  <<
                    " error clonotypes removed)." <<
                    std::endl <<
                    "Parsing is complete. Resulting cloneset size: " <<
                    (size_t) vec.size() <<
                    std::endl;

            return true;
        }


        void parseWordSegment(const string& word, char sep, vector<seg_index_t> &segvec, const GeneSegmentAlphabet& gsa, size_t line_num) {
            stringstream word_stream(word);

            // check if genes are presented in GeneSegmentAlphabet

            if (word_stream.eof()) {
                if (gsa[word].index != 0) {
                    segvec.push_back(gsa[word].index);
                } else {
                    std::cout << "can't find '" << word << "' among genes at line " << line_num << std::endl;
                }
            } else {
                string temp_str = "";
                while (!word_stream.eof()) {
                    getline(word_stream, temp_str, sep);
                    if (gsa[temp_str].index != 0) {
                        segvec.push_back(gsa[temp_str].index);
                    } else {
                        std::cout << "can't find '" << temp_str << "' among genes at line " << line_num << std::endl;
                    }
                }
            }
        }


        void parseWordAlignment(const string& word, char sep, vector<seg_index_t> &segvec, const GeneSegmentAlphabet& gsa, ClonotypeBuilder &clone_builder, char seg) {
            stringstream word_stream(word);
            int aligned_chars = 0;
            if (word_stream.eof()) {
                aligned_chars = stoi(word);
                if (aligned_chars == -1) { aligned_chars = 0; }

                if (seg == 'V') {
                    clone_builder.addValignment(segvec[0], aligned_chars);
                } else {
                    clone_builder.addJalignment(segvec[0], aligned_chars);
                }
            } else {
                string temp_str = "";
                int seg_index = 0;
                while (!word_stream.eof()) {
                    getline(word_stream, temp_str, sep);
                    aligned_chars = stoi(temp_str);
                    if (aligned_chars == -1) { aligned_chars = 0; }

                    if (seg == 'V') {
                        clone_builder.addValignment(segvec[seg_index], aligned_chars);
                    } else {
                        clone_builder.addJalignment(segvec[seg_index], aligned_chars);
                    }

                    ++seg_index;
                }
            }
        }


        string get_prefix(const string &filename) const {
            return "[" + filename + "]: ";
        }

    };
}

#endif