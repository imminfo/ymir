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

#ifndef _GENESEGMENT_H
#define _GENESEGMENT_H


#include <fstream>
#include <sstream>

#include "types.h"
#include "tools.h"


namespace ymir {

    class GeneSegmentAlphabet;
    class VDJRecombinationGenes;


    /**
    * \class GeneSegmentAlphabet
    */
    class GeneSegmentAlphabet {
    public:

        /**
         *
         */
        struct GeneSegment {
            sequence_t allele;
            sequence_t sequence;
            std::string orig_sequence;
            seg_index_t index;

            GeneSegment(const std::string& allele_,
                        const std::string& sequence_,
                        seg_index_t index_)
                    : allele(allele_), sequence(sequence_), orig_sequence(sequence_), index(index_) {
            }
        };


        /**
         *
         */
        GeneSegmentAlphabet(GeneSegments gene_segment, const std::string& name, const std::string& filepath, bool *is_ok = nullptr)
                : _name(name) {
            this->addGeneSegment("other", "", 0);
            bool res = this->read(filepath);
            if (is_ok) {
                *is_ok = res;
            }
            _gene_segment = gene_segment;
        }


        GeneSegmentAlphabet(GeneSegments gene_segment, const std::string& name, const std::vector<std::string>& alleles, const std::vector<std::string>& sequences) {
            this->_name = name;

            this->_vec.reserve(alleles.size());

            this->addGeneSegment("other", "", 0);
            for (size_t i = 0; i < alleles.size(); i++) {
                this->addGeneSegment(alleles[i], sequences[i], i + 1);
            }

            _gene_segment = gene_segment;
        }


        GeneSegmentAlphabet(const GeneSegmentAlphabet &other) {
            _name = other._name;
            _map.clear();
            _map.insert(other._map.begin(), other._map.end());
            _vec.clear();
            _vec.insert(_vec.begin(), other._vec.begin(), other._vec.end());
            _gene_segment = other._gene_segment;
        }


        ~GeneSegmentAlphabet() {}


        const std::string& name() const { return this->_name; }


        const GeneSegment& operator[] (seg_index_t index) const {
            if (index < _vec.size()) { return _vec[index]; }
            return _vec[0];
        }

        const GeneSegment& operator[] (const std::string& name) const {
            if (_map.find(name) != _map.end()) { return _vec[_map.at(name)]; }
            return _vec[0];
        }


        seg_index_t size() const { return this->_vec.size(); }


        seg_index_t max() const { return this->_vec.size() - 1; }


        /**
        * \brief Append P(alindromic)-nucleotides to all gene segment sequences from
        * 3'end, 5'end or both.
        *
        * \param from_5_end Number of P nucleotides to append from the 3' end (i.e., sequences' start).
        * \param from_3_end Number of P nucleotides to append from the 3' end (i.e., sequences' end).
        */
        void appendPalindromicNucleotides(seq_len_t from_5_end, seq_len_t from_3_end) {
            std::string seq_5_end = "", seq_3_end = "";
            for (auto i = 1; i < _vec.size(); ++i) {

                // checks for _vec sequences' length ?

                seq_5_end = "";
                for (auto left_i = 0; left_i < from_5_end; ++ left_i) {
                    seq_5_end = complement(_vec[i].orig_sequence[left_i]) + seq_5_end;
                }

                seq_3_end = "";
                for (auto right_i = 0; right_i < from_3_end; ++right_i) {
                    seq_3_end = seq_3_end + complement(_vec[i].orig_sequence[_vec[i].orig_sequence.size() - 1 - right_i]);
                }

                _vec[i].sequence.reserve(from_5_end + _vec[i].orig_sequence.size() + from_3_end + 2);
                _vec[i].sequence = seq_5_end + _vec[i].orig_sequence + seq_3_end;
            }
        }


        /**
        * \brief Save segment alphabet to the given file as a tab-separated table with 3 columns.
        */
        bool write(const std::string& filepath) const {
            std::ofstream ofs;
            ofs.open(filepath);

            if (ofs.is_open()) {
                ofs << this->_name << '\t' << "Sequences" << std::endl;
                for (seg_index_t i = 1; i < this->_vec.size(); ++i) {
                    ofs << this->_vec[i].allele << '\t' << this->_vec[i].orig_sequence;
                    if (i != this->_vec.size() - 1) {
                        ofs << std::endl;
                    }
                }
            } else {
                return false;
            }

            ofs.close();
            return true;
        }


        bool read(const std::string& filepath) {
            std::ifstream ifs;
            ifs.open(filepath);

            if (ifs.is_open()) {
                std::stringstream line_stream;
                std::string line, gene_name, gene_seq;
                bool header = false;
                seg_index_t cur_index = 1;
                while (!ifs.eof()) {
                    getline(ifs, line);
                    if (line[0] != '\n' && line.size() > 3) {
                        line_stream.str(line);
                        if (header) {
                            // gene segment name
                            getline(line_stream, gene_name, '\t');
                            // gene segment sequence
                            getline(line_stream, gene_seq, '\t');
                            this->addGeneSegment(gene_name, gene_seq, cur_index);
                            ++cur_index;
                        } else {
                            // skip header
                            header = true;
                        }
                        line_stream.clear();
                    }
                }
            } else {
                std::cerr << "Gene segment alphabet [" << this->_name << "] error:" << std::endl << "\tinput file [" << filepath << "] not found" << std::endl;
                return false;
            }

            ifs.close();
            return true;
        }


        seq_len_t maxLength() const {
            seq_len_t res = 0;
            for (auto i = 1; i < _vec.size(); ++i) { 
                res = std::max(res, (seq_len_t) _vec[i].sequence.size());
            }
            return res;
        }


        GeneSegments gene_segment() const { return _gene_segment; }

    protected:

        std::string _name;
        std::unordered_map<std::string, seg_index_t> _map;
        std::vector<GeneSegment> _vec;
        GeneSegments _gene_segment;


        GeneSegmentAlphabet() {}


        void addGeneSegment(const std::string& name, const std::string& seq, seg_index_t index) {
            this->_map[name] = index;
            this->_vec.push_back(GeneSegment(name, seq, index));
        }

    };


    /**
    * \class VDJRecombinationGenes
    */
    class VDJRecombinationGenes {

        typedef unique_ptr<GeneSegmentAlphabet> GSA_ptr;
    public:


        VDJRecombinationGenes() {}


        VDJRecombinationGenes(const std::string& v_segments_name, const std::string& v_segments_file,
                              const std::string& j_segments_name, const std::string& j_segments_file,
                              bool *is_V_ok = nullptr, bool *is_J_ok = nullptr)
            : _V(new GeneSegmentAlphabet(VARIABLE, v_segments_name, v_segments_file, is_V_ok)),
              _J(new GeneSegmentAlphabet(JOINING, j_segments_name, j_segments_file, is_J_ok))
        {
        }


        VDJRecombinationGenes(const std::string& v_segments_name, const std::string& v_segments_file,
                              const std::string& j_segments_name, const std::string& j_segments_file,
                              const std::string& d_segments_name, const std::string& d_segments_file,
                              bool *is_V_ok = nullptr, bool *is_J_ok = nullptr, bool *is_D_ok = nullptr)
            : _V(new GeneSegmentAlphabet(VARIABLE, v_segments_name, v_segments_file, is_V_ok)),
              _J(new GeneSegmentAlphabet(JOINING, j_segments_name, j_segments_file, is_J_ok)),
              _D(new GeneSegmentAlphabet(DIVERSITY, d_segments_name, d_segments_file, is_D_ok))
        {
        }


        VDJRecombinationGenes(const std::string& V_name, const std::vector<std::string>& V_alleles, const std::vector<std::string>& V_sequences,
                              const std::string& J_name, const std::vector<std::string>& J_alleles, const std::vector<std::string>& J_sequences)
            : _V(new GeneSegmentAlphabet(VARIABLE, V_name, V_alleles, V_sequences)),
              _J(new GeneSegmentAlphabet(JOINING, J_name, J_alleles, J_sequences))
        {
        }


        VDJRecombinationGenes(const std::string& V_name, const std::vector<std::string>& V_alleles, const std::vector<std::string>& V_sequences,
                              const std::string& J_name, const std::vector<std::string>& J_alleles, const std::vector<std::string>& J_sequences,
                              const std::string& D_name, const std::vector<std::string>& D_alleles, const std::vector<std::string>& D_sequences)
            : _V(new GeneSegmentAlphabet(VARIABLE, V_name, V_alleles, V_sequences)),
              _J(new GeneSegmentAlphabet(JOINING, J_name, J_alleles, J_sequences)),
              _D(new GeneSegmentAlphabet(DIVERSITY, D_name, D_alleles, D_sequences))
        {
        }


        VDJRecombinationGenes(const VDJRecombinationGenes &other) {
            _V.reset(new GeneSegmentAlphabet(*other._V));
            _J.reset(new GeneSegmentAlphabet(*other._J));
            if (other._D) {
                _D.reset(new GeneSegmentAlphabet(*other._D));
            } else {
                _D.release();
            }
        }


        VDJRecombinationGenes& operator=(const VDJRecombinationGenes &other) {
            _V.reset(new GeneSegmentAlphabet(*other._V));
            _J.reset(new GeneSegmentAlphabet(*other._J));
            if (other._D) {
                _D.reset(new GeneSegmentAlphabet(*other._D));
            } else {
                _D.release();
            }
        }


        ~VDJRecombinationGenes()
        {
        }


        bool is_vdj() const { return this->_D != nullptr; }


        /**
         *
         */
        ///@{
        const GeneSegmentAlphabet& V() const { return *(this->_V); }

        const GeneSegmentAlphabet& J() const { return *(this->_J); }

        const GeneSegmentAlphabet& D() const { return *(this->_D); }
        ///@}


        /**
        * \param segindex 0 for V(ariable), 1 for J(oining), 2 for D(iversity) gene segment.
        */
        void appendPalindromicNucleotides(GeneSegments gene, seq_len_t from_5_end = 4, seq_len_t from_3_end = 4) {
            if (gene == VARIABLE) {
                if (this->_V) { this->_V->appendPalindromicNucleotides(from_5_end, from_3_end); }
            } else if (gene == JOINING) {
                if (this->_J) { this->_J->appendPalindromicNucleotides(from_5_end, from_3_end); }
            } else if (gene == DIVERSITY) {
                if (this->_D) { this->_D->appendPalindromicNucleotides(from_5_end, from_3_end); }
            }
        }


        bool write(const std::string& v_segments_file = "vsegments.txt",
                const std::string& j_segments_file = "jsegments.txt",
                const std::string& d_segments_file = "dsegments.txt") const {
            if (this->_V->write(v_segments_file)) {
                if (this->_J->write(j_segments_file)) {
                    if (!this->is_vdj()) {
                        return true;
                    } else {
                        return this->_D->write(d_segments_file);
                    }
                }
            }
            return false;
        }

    protected:

        GSA_ptr _V, _J, _D;

    };
}

#endif