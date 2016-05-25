//
// Created by Vadim N. on 09/05/2016.
//

#ifndef YMIR_PARSER_BASE_H
#define YMIR_PARSER_BASE_H


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
            OVERWRITE,        // Try to align all gene segments
            USE_PROVIDED,     // Use information about aligned gene segments - start / end / length (no errors)
            REALIGN_PROVIDED  // Re-align aligned gene segments using start / end (with errors)
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
    template <typename ClonotypeType>
    class ParserBase {

    public:

        typedef ClonotypeType clonotype_t;

        typedef VDJAlignerBase<clonotype_t> vdj_aligner_t;

        typedef ClonotypeVector<clonotype_t> clonotype_vector_t;

        typedef Cloneset<clonotype_t> cloneset_t;

        /**
        * \typedef ParserConfig
        *
        * \brief Parameters of this parser: parser's name, names of the columns, separators, sequences and other.
        */
        typedef Json::Value ParserConfig;


        ParserBase() {
//            _config_is_loaded = false;
        }


        ParserBase(vdj_aligner_t *aligner)
                : _aligner(aligner)
        {
        }


        ~ParserBase()
        {
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

            if (!_aligner) {
                std::cout << "Repertoire parser error:" << "\tsupply the aligner to the parser." << endl;
                return false;
            }

            _stream.open(filepath);
            if (_stream.good()) {
                std::cout << "Open [" << filepath << "] for reading" << endl;
                _genes = gene_segments;
                std::cout << "-- gene segments assigned to the parser" << endl;
                _recomb = recomb;
                _opts = opts;
                std::cout << "-- options assigned" << endl;
                _aligner->set_genes(_genes);
                std::cout << "-- gene segments assigned to the aligner" << endl;
                _aligner->set_parameters(params);
                std::cout << "-- parameters assigned to the aligner" << endl;
                _status = true;
                return true;
            } else {
                std::cout << "Repertoire parser error:" << "\tinput file [" << filepath << "] not found" << endl;
            }

            return false;
        }


        /**
         * \param cloneset Pointer to clonal repertoire object to which data will be uploaded.
         */
        bool parse(cloneset_t *cloneset, size_t max_clonotype_count = (size_t)-1) {
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

            clonotype_vector_t clonevec;
            clonevec.reserve(DEFAULT_REPERTOIRE_RESERVE_SIZE);

            this->parseRepertoire(clonevec, max_clonotype_count);
            cloneset->swap(clonevec);

            return true;
        }


        bool openAndParse(const std::string &filepath,
                          cloneset_t *cloneset,
                          const VDJRecombinationGenes &gene_segments,
                          Recombination recomb,
                          AlignmentColumnOptions opts = AlignmentColumnOptions().setV(AlignmentColumnOptions::USE_PROVIDED).setJ(AlignmentColumnOptions::USE_PROVIDED).setD(AlignmentColumnOptions::OVERWRITE),
                          VDJAlignerParameters params = VDJAlignerParameters()) {
            if (this->open(filepath, gene_segments, recomb, opts, params)) {
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
        RepertoireParserStatistics _stats;
        unique_ptr<vdj_aligner_t> _aligner;
        bool _status;
        bool _read_header;


        void parseRepertoire(clonotype_vector_t& vec, size_t max_clonotype_count)
        {
            char column_sep ='\t',
                    segment_sep = ',',
                    internal_sep = '|',
                    alignment_sep = ';',
                    start_bracket = '(',
                    end_bracket = ')';

            AlignmentColumnOptions::Action align_V_opt = _opts.align_V,
                    align_J_opt = _opts.align_J,
                    align_D_opt = _opts.align_D;

            stringstream column_stream, symbol_stream, temp_stream;
            string line, segment_word, sequence;

            int clonotype_num = 0, line_num = _stats.count_all + 1;

            bool align_ok = true;

            vector<seg_index_t> vseg, jseg, dseg;
            string temp_str;

            _aligner->setRecombination(_recomb);

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
                    if (clonotype_t::sequence_type == NUCLEOTIDE) {
                        getline(column_stream, sequence, column_sep);
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                        getline(column_stream, sequence, column_sep);
                    }
                    _aligner->setSequence(sequence);

                    //
                    // Parse Variable genes
                    //
                    if (align_V_opt == OVERWRITE) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseSegment(symbol_stream, segment_word, vseg, VARIABLE, _genes.V(), line_num, segment_sep, temp_str);
                    }

                    //
                    // Parse Diversity genes
                    //
                    if (_recomb == VJ_RECOMB) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        if (align_D_opt == OVERWRITE) {
                            column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                        } else {
                            getline(column_stream, segment_word, column_sep);
                            this->parseSegment(symbol_stream, segment_word, dseg, DIVERSITY, _genes.D(), line_num, segment_sep, temp_str);
                        }
                    }

                    //
                    // Parse Joining genes
                    //
                    if (align_J_opt == OVERWRITE) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                    } else {
                        getline(column_stream, segment_word, column_sep);
                        this->parseSegment(symbol_stream, segment_word, jseg, JOINING, _genes.J(), line_num, segment_sep, temp_str);
                    }

                    //
                    // Parse Variable alignments
                    //
                    if (do_align_V) {
                        column_stream.ignore(numeric_limits<streamsize>::max(), column_sep);
                        if (!_aligner->alignVar()) {
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
                            if (!_aligner->alignDiv()) {
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
                        if (!_aligner->alignJoi()) {
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
                    vec.push_back(_aligner->buildClonotype());
                    std::cout << vec[vec.size() - 1].toString(_genes) << std::endl;
                }
            }
        }


        virtual void parseSegment(stringstream &symbol_stream,
                                  const string &segment_word,
                                  vector<seg_index_t> &segvec,
                                  GeneSegments gene,
                                  const GeneSegmentAlphabet &gsa,
                                  size_t line_num,
                                  char segment_sep,
                                  string &temp_str) = 0;


        virtual void parseAlignment(stringstream &symbol_stream,
                            const string &segment_word,
                            const vector<seg_index_t> &segvec,
                            GeneSegments gene,
                            const GeneSegmentAlphabet &gsa,
                            size_t line_num,
                            char segment_sep,
                            char internal_sep,
                            string &temp_str,
                            stringstream &temp_stream) = 0;


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

}

#endif //YMIR_PARSER_BASE_H
