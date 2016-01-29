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
#include "errcorr.h"
#include "repertoire.h"
#include "genesegment.h"


namespace ymir {

    /**
    * \class RepertoireParser
    *
    * \brief Parser for text files with repertoire data. By default it's MiTCR parser.
    * To make new parsers inherit from this class and rewrite virtual private method
    * "parseRepertoire".
    */
    class RepertoireParser {

    public:

        static const size_t default_block_size = 100000;


        /**
        * \enum ALIGNMENT_COLUMN_ACTION
        *
        * \brief Specify an action to perform with an alignment column:
        * either overwrite found alignments (OVERWRITE)
        * or use found (USE_PROVIDED) alignments in the input file.
        */
        enum ALIGNMENT_COLUMN_ACTION {
            OVERWRITE,
            USE_PROVIDED
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
        bool parse(const std::string& filepath,
                   Cloneset *rep,
                   const VDJRecombinationGenes& gene_segments,
                   SequenceType seq_type,
                   Recombination recomb,
                   AlignmentColumnOptions opts = AlignmentColumnOptions().setV(USE_PROVIDED).setJ(USE_PROVIDED).setD(OVERWRITE),
                   const AbstractAligner& aligner = NaiveNucleotideAligner()) {

            if (recomb == UNDEF_RECOMB) {
                std::cout << "Repertoire parser error:" << "\tno recombination type for [" << filepath << "]" << endl;
                return false;
            }

            ClonotypeVector clonevec;
            clonevec.reserve(DEFAULT_REPERTOIRE_RESERVE_SIZE);

            ifstream ifs;
            ifs.open(filepath);
            bool res = false;
            if (ifs.is_open()) {
                std::cout << "Parsing input file:\t" << filepath << endl;
                res = this->parseRepertoire(filepath, 
                                            ifs, 
                                            clonevec, 
                                            gene_segments, 
                                            aligner, 
                                            seq_type, 
                                            opts, 
                                            recomb);
                if (res) { 
                    rep->swap(clonevec);
                }
            } else {
                std::cout << "Repertoire parser error:" << "\tinput file [" << filepath << "] not found" << endl;
                res = false;
            }

            return res;
        }


        /**
         * \brief Open the parser in the stream mode - iteratively parse an
         * input repertoire by blocks and return each block.
         */
        bool stream(const std::string& filepath, 
                    const VDJRecombinationGenes& gene_segments,
                    SequenceType seq_type,
                    Recombination recomb,
                    AlignmentColumnOptions opts = AlignmentColumnOptions().setV(USE_PROVIDED).setJ(USE_PROVIDED).setD(OVERWRITE),
                    const AbstractAligner& aligner = NaiveNucleotideAligner())
        {
            if (recomb == UNDEF_RECOMB) {
                std::cout << "Repertoire parser error:" << "\tno recombination type for [" << filepath << "]" << endl;
                return false;
            }

            _stream.open(filepath);
            if (_stream.is_open()) {
                std::cout << "Open the stream to the input file:\t" << filepath << endl;
                return true;
            } else {
                std::cout << "Repertoire parser error:" << "\tinput file [" << filepath << "] not found" << endl;
                return false;
            }
        }


        /**
         * \brief Get the next Block from the parser stream.
         */
        void nextBlock(Cloneset *rep, size_t block_size = default_block_size) {
            if (_stream.good()) {
                ClonotypeVector clonevec;
                clonevec.reserve(DEFAULT_REPERTOIRE_RESERVE_SIZE);

                bool res = this->parseRepertoire(filepath, 
                                                ifs, 
                                                clonevec, 
                                                gene_segments, 
                                                aligner, 
                                                seq_type, 
                                                opts, 
                                                recomb, 
                                                block_size);
                if (res) {
                    rep->swap(clonevec);
                }
            } else {
                std::cout << "Repertoire parser error: bad / closed stream to [" << filepath << "]" << endl;
            }

            if (_stream.eof()) { 
                _stream.close()
            }
        }


    protected:

//        ParserConfig _config;
//        bool _config_is_loaded;
        std::ifstream _stream;


        virtual bool parseRepertoire(const string &filename,
                                     ifstream& ifs,
                                     ClonotypeVector& vec,
                                     const VDJRecombinationGenes& gene_segments,
                                     const AbstractAligner& aligner,
                                     SequenceType seq_type,
                                     AlignmentColumnOptions opts,
                                     Recombination recomb, 
                                     size_t clonotype_count = (size_t)-1)
        {
            char column_sep ='\t',
                 segment_sep = ',',
                 internal_sep = '|',
                 alignment_sep = ';',
                 start_bracket = '(',
                 end_bracket = ')';

            bool do_align_V = opts.align_V == OVERWRITE,
                 do_align_J = opts.align_J == OVERWRITE,
                 do_align_D = opts.align_D == OVERWRITE;

            stringstream column_stream, symbol_stream, temp_stream;
            string line, segment_word, sequence;

            int index = 0;

            vector<seg_index_t> vseg, jseg, dseg;
            string temp_str;

            ClonotypeBuilder clone_builder;

            int glob_index = 1, bad_index = 0;

            // Skip header
            getline(ifs, line);
            while (!ifs.eof() && glob_index < clonotype_count) {
                // Start processing clonotypes
                getline(ifs, line);
                if (line.size() > 2) {
                    // parse body and build clonotypes from each line
                    vseg.clear();
                    jseg.clear();
                    dseg.clear();

                    column_stream.clear();
                    column_stream.str(line);

                    clone_builder.setRecombination(recomb);

                    // Get nucleotide or amino acid sequence
                    if (seq_type == NUCLEOTIDE) {
                        getline(column_stream, sequence, column_sep);
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                        getline(column_stream, sequence, column_sep);
                    }

                    clone_builder.setSequence(sequence);
                    clone_builder.setSequenceType(seq_type);


                    // Parse Variable genes
                    if (do_align_V) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseSegment(symbol_stream, segment_word, vseg, gene_segments.V(), clone_builder, glob_index, segment_sep, temp_str);
                    }

                    // Parse Diversity genes
                    if (recomb == VJ_RECOMB) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        if (do_align_D) {
                            column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                        } else {
                            getline(column_stream, segment_word, column_sep);
                            this->parseSegment(symbol_stream, segment_word, dseg, gene_segments.D(), clone_builder, glob_index, segment_sep, temp_str);
                        }
                    }

                    // Parse Joining genes
                    if (do_align_J) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseSegment(symbol_stream, segment_word, jseg, gene_segments.J(), clone_builder, glob_index, segment_sep, temp_str);
                    }


                    // Parse Variable alignments
                    if (do_align_V) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);

                        //
                        // alignment here
                        //
                        // TODO: implement V alignment in parser
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseAlignment(symbol_stream, segment_word, vseg, gene_segments.V(), clone_builder, glob_index, segment_sep, internal_sep, temp_str, temp_stream);
                    }

                    bool align_ok = true;
                    // Parse Diversity alignments
                    if (recomb == VJ_RECOMB) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        if (do_align_D) {
                            column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);

                            //
                            // alignment here
                            //
                            // TODO: implement D alignment in parser
                            align_ok = false;
                            for (seg_index_t seg_i = 1; seg_i <= gene_segments.D().max(); ++seg_i) {
                                AbstractAligner::LocalAlignmentIndices indices =
                                        aligner.alignLocal(gene_segments.D()[seg_i].sequence,
                                                           sequence,
                                                           DEFAULT_DIV_GENE_MIN_LEN);
                                if (indices.size()) {
                                    align_ok = true;
                                    for (size_t align_i = 0; align_i < indices.size(); ++align_i) {
                                        clone_builder.addDivAlignment(seg_i, indices[align_i]);
                                    }
                                }
                            }
                            if (!align_ok) {
                                cout << "Diversity gene could NOT be aligned with the given minimal gene length (min gene length " << (size_t) DEFAULT_DIV_GENE_MIN_LEN << ", line " << (size_t) glob_index << ")" << endl;
                                ++bad_index;
                            }
                        } else {
                            getline(column_stream, segment_word, column_sep);
                            this->parseAlignment(symbol_stream, segment_word, dseg, gene_segments.D(), clone_builder, glob_index, segment_sep, internal_sep, temp_str, temp_stream);
                        }
                    }

                    // Parse Joining alignments
                    if (do_align_J) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);

                        //
                        // alignment here
                        //
                        // TODO: implement J alignment in parser
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseAlignment(symbol_stream, segment_word, jseg, gene_segments.J(), clone_builder, glob_index, segment_sep, internal_sep, temp_str, temp_stream);
                    }

                    ++index;

                    if (glob_index % 50000 == 0) {
                        cout << this->get_prefix(filename) + "parsed " << (size_t) glob_index << " lines" << endl;
                    }
                    ++glob_index;

                    //
                    // remove bad clonotypes here ???
                    //
                    if (align_ok) {
                        vec.push_back(clone_builder.buildClonotype());
                    }

//                    vec.push_back(clone_builder.buildClonotype());
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


        void parseSegment(stringstream &symbol_stream,
                          const string &segment_word,
                          vector<seg_index_t> &segvec,
                          const GeneSegmentAlphabet &gsa,
                          ClonotypeBuilder &clone_builder,
                          size_t line_num,
                          char segment_sep,
                          string &temp_str)
        {
            symbol_stream.clear();
            symbol_stream.str(segment_word);

            if (symbol_stream.eof()) {
                if (gsa[segment_word].index != 0) {
                    segvec.push_back(gsa[segment_word].index);
                } else {
                    std::cout << "can't find '" << segment_word << "' among gene segments at the line " << (size_t) line_num << std::endl;
                }
            } else {
                while (!symbol_stream.eof()) {
                    getline(symbol_stream, temp_str, segment_sep);
                    if (gsa[temp_str].index != 0) {
                        segvec.push_back(gsa[temp_str].index);
                    } else {
                        std::cout << "can't find '" << temp_str << "' among genes at line " << line_num << std::endl;
                    }
                }
            }
        }


        void parseAlignment(stringstream &symbol_stream,
                            const string &segment_word,
                            const vector<seg_index_t> &segvec,
                            const GeneSegmentAlphabet &gsa,
                            ClonotypeBuilder &clone_builder,
                            size_t line_num,
                            char segment_sep,
                            char internal_sep,
                            string &temp_str,
                            stringstream &temp_stream)
        {
            symbol_stream.clear();
            symbol_stream.str(segment_word);

            int gene_start, seq_start, alignment_len;
            seg_index_t seg_order = 0;

            if (symbol_stream.eof()) {
                getline(symbol_stream, temp_str, segment_sep);

                temp_stream.clear();
                temp_stream.str(temp_str);

                getline(temp_stream, temp_str, internal_sep);
                gene_start = std::atoi(temp_str.c_str());
                getline(temp_stream, temp_str, internal_sep);
                seq_start = std::atoi(temp_str.c_str());
                getline(temp_stream, temp_str, internal_sep);
                alignment_len = std::atoi(temp_str.c_str());

                clone_builder.addAlignment(gsa.gene_segment(), segvec[0], gene_start, seq_start, alignment_len);

            } else {
                while (!symbol_stream.eof()) {
                    getline(symbol_stream, temp_str, segment_sep);

                    temp_stream.clear();
                    temp_stream.str(temp_str);

                    getline(temp_stream, temp_str, internal_sep);
                    gene_start = std::atoi(temp_str.c_str());

                    getline(temp_stream, temp_str, internal_sep);
                    seq_start = std::atoi(temp_str.c_str());

                    getline(temp_stream, temp_str, internal_sep);
                    alignment_len = std::atoi(temp_str.c_str());
                    // "0" in gene_start means that there is no information from what letter
                    // the J gene segment was started, so Ymir by default will compute it
                    // assuming that J segment is aligned at the very end of the input sequence.
                    if (gene_start == 0) {
                        gene_start = gsa[segvec[seg_order]].sequence.size() - alignment_len + 1;
                    }

                    clone_builder.addAlignment(gsa.gene_segment(), segvec[seg_order], gene_start, seq_start, alignment_len);

                    ++seg_order;
                }
            }
        }


        std::string get_prefix(const string &filename) const {
            return "[" + filename + "]: ";
        }

    };
}

#endif