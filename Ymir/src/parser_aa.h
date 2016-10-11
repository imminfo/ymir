//
// Created by Vadim N. on 09/05/2016.
//

#ifndef YMIR_PARSER_AA_H
#define YMIR_PARSER_AA_H

#include "parser_base.h"


namespace ymir {

    class ParserAA : public ParserBase<ClonotypeAA> {

    public:

        ParserAA() {
//            _config_is_loaded = false;
        }


        ParserAA(vdj_aligner_t *aligner) : ParserBase(aligner)
        {
        }

    protected:

        void parseSegment(stringstream &symbol_stream,
                          const string &segment_word,
                          vector<seg_index_t> &segvec,
                          GeneSegments gene,
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
                    switch (gene) {
                        case VARIABLE: {
                            ++_stats.bad_V_seg;
                            break;
                        }
                        case DIVERSITY: {
                            ++_stats.bad_D_seg;
                            break;
                        }
                        case JOINING: {
                            ++_stats.bad_J_seg;
                            break;
                        }
                        default: { }
                    }
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

                _aligner->addAlignment(gene, segvec[seg_order], gene_start, seq_start, alignment_len);

                ++seg_order;
            }
        }

    };

}

#endif //YMIR_PARSER_AA_H
