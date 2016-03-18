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
#include "clonotype_builder.h"
#include "errcorr.h"
#include "repertoire.h"


using std::getline;


namespace ymir {

    /**
     *
     */
    struct RepertoireParserStatistics {

        RepertoireParserStatistics()
        {
            this->reset();
        }

        void reset() {
            count_all = 0;
            bad_V_seg = 0;
            bad_D_seg = 0;
            bad_J_seg = 0;
            no_V_algn = 0;
            no_D_algn = 0;
            no_J_algn = 0;
            bad_V_len = 0;
            bad_J_len = 0;
        }

        void print() {
            std::cout <<
                    "Parsing complete. Parsed " << (size_t) count_all << " lines, found clonotypes with:" << std::endl <<
                    "\tunrecognised V segments:\t" << (size_t) bad_V_seg << std::endl <<
                    "\tunrecognised D segments:\t" << (size_t) bad_D_seg << std::endl <<
                    "\tunrecognised J segments:\t" << (size_t) bad_J_seg << std::endl <<
                    "\tunaligned V segments:\t" << (size_t) no_V_algn << std::endl <<
                    "\tunaligned D segments:\t" << (size_t) no_D_algn << std::endl <<
                    "\tunaligned J segments:\t" << (size_t) no_J_algn << std::endl <<
                    "\tbad alignment length of V segments:\t" << (size_t) bad_V_len << " (repaired)" << std::endl <<
                    "\tbad alignment length of J segments:\t" << (size_t) bad_J_len << " (repaired)" << std::endl <<
                    "Resulting cloneset size: " << (size_t) count_all << std::endl;
        }


        template <GeneSegments GENE>
        void update_bad_seg();


        template <GeneSegments GENE>
        void update_no_algn();


        template <GeneSegments GENE>
        void update_bad_len();


        size_t count_all, bad_V_seg, bad_D_seg, bad_J_seg, no_V_algn, no_D_algn, no_J_algn, bad_V_len, bad_J_len;

    };


//        enum AlignmentColumnsAction {
//            USE_PROVIDED,
//            ALIGN_ONLY_D,
//            ALIGN_ALL
//        };


    struct AlignmentColumnOptions {

        /**
         * \enum ALIGNMENT_COLUMN_ACTION
         *
         * \brief Specify an action to perform with an alignment column:
         * either overwrite found alignments (OVERWRITE)
         * or use found (USE_PROVIDED) alignments in the input file.
        */
        enum Action {
            OVERWRITE,
            USE_PROVIDED
        };

        Action align_V, align_J, align_D;

        AlignmentColumnOptions() {}


        AlignmentColumnOptions(Action V, Action J)
            : align_V(V), align_J(J)
        {
        }


        AlignmentColumnOptions(Action V, Action D, Action J)
            : align_V(V), align_D(D), align_J(J)
        {
        }


        AlignmentColumnOptions& setV(Action action) {
            this->align_V = action;
            return *this;
        }

        AlignmentColumnOptions& setJ(Action action) {
            this->align_J = action;
            return *this;
        }

        AlignmentColumnOptions& setD(Action action) {
            this->align_D = action;
            return *this;
        }
    };


    /**
    * \class RepertoireParser
    *
    * \brief Parser for text files with repertoire data. By default it's MiTCR parser.
    * To make new parsers inherit from this class and rewrite virtual private method
    * "parseRepertoire".
    */
    template <typename Aligner >
    class RepertoireParser {

    public:

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
        * \param filepath Path to a file with sequences.
        * \param gene_segments Alphabets of gene segments.
        * \param seq_type
        * \param recomb
        * \param opts What action to do with columns with V, D and J alignments.
        *
        * \return True if found has been successfully processed, false otherwise.
        */
        bool open(const std::string &filepath,
                  const VDJRecombinationGenes &gene_segments,
                  SequenceType seq_type,
                  Recombination recomb,
                  AlignmentColumnOptions opts = AlignmentColumnOptions().setV(AlignmentColumnOptions::USE_PROVIDED).setJ(AlignmentColumnOptions::USE_PROVIDED).setD(AlignmentColumnOptions::OVERWRITE),
                  VDJAlignerParameters params = VDJAlignerParameters()) {
            if (_stream.is_open()) { _stream.close(); }

            _status = false;
            _read_header = true;
            _stats.reset();

            if (recomb == UNDEF_RECOMB) {
                std::cout << "Repertoire parser error:" << "\tno recombination type for [" << filepath << "]" << endl;
                return false;
            }

            _stream.open(filepath);
            if (_stream.good()) {
                std::cout << "Open [" << filepath << "] for reading" << endl;
                _genes = gene_segments;
                std::cout << "-- gene segments assigned to the parser" << endl;
                _seq_type = seq_type;
                _recomb = recomb;
                _opts = opts;
                std::cout << "-- options assigned" << endl;
                _status = true;
                _aligner.set_genes(_genes);
                std::cout << "-- gene segments assigned to the aligner" << endl;
                _aligner.set_parameters(params);
                std::cout << "-- parameters assigned to the aligner" << endl;
                return true;
            } else {
                std::cout << "Repertoire parser error:" << "\tinput file [" << filepath << "] not found" << endl;
            }

            return false;
        }


        /**
         * \param cloneset Pointer to clonal repertoire object to which data will be uploaded.
         */
        bool parse(Cloneset *cloneset, size_t max_clonotype_count = (size_t)-1) {
            if (_stream.eof()) {
                _stats.print();
                _stats.reset();
                _status = false;
                _stream.close();
                return false;
            }

            if (!_status) {
                std::cout << "Something bad is happening - can\'t parse the input file. Perhaps you need to open it first with open()?" << endl;
                return false;
            }

            ClonotypeVector clonevec;
            clonevec.reserve(DEFAULT_REPERTOIRE_RESERVE_SIZE);

            this->parseRepertoire(clonevec, max_clonotype_count);
            cloneset->swap(clonevec);

            return true;
        }


        bool openAndParse(const std::string &filepath,
                          Cloneset *cloneset,
                          const VDJRecombinationGenes &gene_segments,
                          SequenceType seq_type,
                          Recombination recomb,
                          AlignmentColumnOptions opts = AlignmentColumnOptions().setV(AlignmentColumnOptions::USE_PROVIDED).setJ(AlignmentColumnOptions::USE_PROVIDED).setD(AlignmentColumnOptions::OVERWRITE),
                          VDJAlignerParameters params = VDJAlignerParameters()) {
            if (this->open(filepath, gene_segments, seq_type, recomb, opts, params)) {
                this->parse(cloneset);
                _stats.print();
                if (_stream.is_open()) { _stream.close(); }
                return true;
            }

            return false;
        }

    protected:

        std::ifstream _stream;
        VDJRecombinationGenes _genes;
        Recombination _recomb;
        AlignmentColumnOptions _opts;
        SequenceType _seq_type;
        RepertoireParserStatistics _stats;
        Aligner _aligner;
        bool _status;
        bool _read_header;


        void parseRepertoire(ClonotypeVector& vec, size_t max_clonotype_count)
        {
            char column_sep ='\t',
                 segment_sep = ',',
                 internal_sep = '|',
                 alignment_sep = ';',
                 start_bracket = '(',
                 end_bracket = ')';

            bool do_align_V = _opts.align_V == AlignmentColumnOptions::OVERWRITE,
                 do_align_J = _opts.align_J == AlignmentColumnOptions::OVERWRITE,
                 do_align_D = _opts.align_D == AlignmentColumnOptions::OVERWRITE;

            stringstream column_stream, symbol_stream, temp_stream;
            string line, segment_word, sequence;

            int clonotype_num = 0, line_num = _stats.count_all + 1;

            bool align_ok = true;

            vector<seg_index_t> vseg, jseg, dseg;
            string temp_str;

            _aligner.setSequenceType(_seq_type);
            _aligner.setRecombination(_recomb);

            // Skip header
            if (_read_header) {
                getline(_stream, line);
                _read_header = false;
            }
            while (!_stream.eof() && clonotype_num < max_clonotype_count) {
                // Start processing clonotypes
                getline(_stream, line);

                if (line.size() > 2) {
                    // parse body and build clonotypes from each line
                    vseg.clear();
                    jseg.clear();
                    dseg.clear();

                    column_stream.clear();
                    column_stream.str(line);

                    //
                    // Get nucleotide or amino acid sequence
                    //
                    if (_seq_type == NUCLEOTIDE) {
                        getline(column_stream, sequence, column_sep);
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                        getline(column_stream, sequence, column_sep);
                    }
                    _aligner.setSequence(sequence);

                    //
                    // Parse Variable genes
                    //
                    if (do_align_V) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseSegment<VARIABLE>(symbol_stream, segment_word, vseg, _genes.V(), line_num, segment_sep, temp_str);
                    }

                    //
                    // Parse Diversity genes
                    //
                    if (_recomb == VJ_RECOMB) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        if (do_align_D) {
                            column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                        } else {
                            getline(column_stream, segment_word, column_sep);
                            this->parseSegment<DIVERSITY>(symbol_stream, segment_word, dseg, _genes.D(), line_num, segment_sep, temp_str);
                        }
                    }

                    //
                    // Parse Joining genes
                    //
                    if (do_align_J) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseSegment<JOINING>(symbol_stream, segment_word, jseg, _genes.J(), line_num, segment_sep, temp_str);
                    }

                    //
                    // Parse Variable alignments
                    //
                    if (do_align_V) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                        if (!_aligner.alignVar()) {
                            _stats.update_no_algn<VARIABLE>();
                        }
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseAlignment(symbol_stream, segment_word, vseg, VARIABLE, _genes.V(), line_num, segment_sep, internal_sep, temp_str, temp_stream);
                    }

                    //
                    // Parse Diversity alignments
                    //
                    if (_recomb == VJ_RECOMB) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        if (do_align_D) {
                            column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                            if (!_aligner.alignDiv()) {
                                _stats.update_no_algn<DIVERSITY>();
                            }
                        } else {
                            getline(column_stream, segment_word, column_sep);
                            this->parseAlignment(symbol_stream, segment_word, dseg, DIVERSITY, _genes.D(), line_num, segment_sep, internal_sep, temp_str, temp_stream);
                        }
                    }

                    //
                    // Parse Joining alignments
                    //
                    if (do_align_J) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                        if (!_aligner.alignJoi()) {
                            _stats.update_no_algn<JOINING>();
                        }
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseAlignment(symbol_stream, segment_word, jseg, JOINING, _genes.J(), line_num, segment_sep, internal_sep, temp_str, temp_stream);
                    }

                    ++clonotype_num;
                    ++_stats.count_all;
                    if (_stats.count_all % 50000 == 0) {
                        cout << "Parsed " << (size_t) _stats.count_all << " lines" << endl;
                    }

                    //
                    // TODO: remove bad clonotypes here ???
                    //
                    vec.push_back(_aligner.buildClonotype());
                }
            }
        }


        template <GeneSegments GENE>
        void parseSegment(stringstream &symbol_stream,
                          const string &segment_word,
                          vector<seg_index_t> &segvec,
                          const GeneSegmentAlphabet &gsa,
                          size_t line_num,
                          char segment_sep,
                          string &temp_str)
        {
            symbol_stream.clear();
            symbol_stream.str(segment_word);

            while (!symbol_stream.eof()) {
                getline(symbol_stream, temp_str, segment_sep);
                if (gsa[temp_str].index != 0) {
                    segvec.push_back(gsa[temp_str].index);
                } else {
                    _stats.update_bad_seg<GENE>();
                }
            }
        }


        void parseAlignment(stringstream &symbol_stream,
                            const string &segment_word,
                            const vector<seg_index_t> &segvec,
                            GeneSegments gene,
                            const GeneSegmentAlphabet &gsa,
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

                if (alignment_len > gsa[segvec[seg_order]].sequence.size()) {
                    alignment_len = gsa[segvec[seg_order]].sequence.size();
                    if (gene == VARIABLE) {
                        ++_stats.bad_V_len;
                    } else if (gene == JOINING) {
                        ++_stats.bad_J_len;
                    }
                }

                // "0" in gene_start means that there is no information from what letter
                // the J gene segment was started, so Ymir by default will compute it
                // assuming that J segment is aligned at the very end of the input sequence.
                if (gene_start == 0) {
                    gene_start = gsa[segvec[seg_order]].sequence.size() - alignment_len + 1;
                }

                _aligner.addAlignment(gene, segvec[seg_order], gene_start, seq_start, alignment_len);

                ++seg_order;
            }
        }


        std::string get_prefix(const string &filename) const {
            return "[" + filename + "]: ";
        }

    };


    template <>
    void RepertoireParserStatistics::update_bad_seg<VARIABLE>() { ++bad_V_seg; }


    template <>
    void RepertoireParserStatistics::update_bad_seg<DIVERSITY>() { ++bad_D_seg; }


    template <>
    void RepertoireParserStatistics::update_bad_seg<JOINING>() { ++bad_J_seg; }


    template <>
    void RepertoireParserStatistics::update_no_algn<VARIABLE>() { ++no_V_algn; }


    template <>
    void RepertoireParserStatistics::update_no_algn<DIVERSITY>() { ++no_D_algn; }


    template <>
    void RepertoireParserStatistics::update_no_algn<JOINING>() { ++no_J_algn; }

//
//    template <GeneSegments GENE>
//    inline void add_alignment(Aligner &aligner, seg_index_t seg_index, seq_len_t genestart, seq_len_t seqstart, seq_len_t alignment_len);
//
//
//    inline void add_alignment<Aligner, VARIABLE>(Aligner &aligner, seg_index_t seg_index, seq_len_t genestart, seq_len_t seqstart, seq_len_t alignment_len) {
//        aligner.addVarAlignment(seg_index, genestart, seqstart, alignment_len);
//    };



    typedef RepertoireParser<NaiveCDR3NucleotideAligner> NaiveNucParser;


    typedef RepertoireParser<CDR3NucleotideAligner> CDR3NucParser;


    typedef RepertoireParser<SmithWatermanNoGapAligner> SWNGParser;


    typedef RepertoireParser<SmithWatermanAligner> SWParser;


//    typedef RepertoireFileParser<NoGapAlignmentVector> RepertoireParser;
//
//    typedef RepertoireFileParser<GappedAlignmentVector> ErrorCorrectingParser;

}

#endif