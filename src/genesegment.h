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


#include "iostream"
#include "fstream"
#include "sstream"
#include "unordered_map"
#include "vector"

#include "types.h"


using namespace std;


namespace ymir {

    class GeneSegmentAlphabet;
    class VDJRecombinationGenes;


    /**
    * \class GeneSegmentAlphabet
    */
    class GeneSegmentAlphabet {
    public:

        struct GeneSegment {
            string allele;
            string sequence;
            segindex_t index;

            GeneSegment(const string& allele_, const string& sequence_, segindex_t index_)
                    : allele(allele_), sequence(sequence_), index(index_) {
            }
        };


        GeneSegmentAlphabet(const string& name, const string& filepath, bool *is_ok = nullptr)
                : _name(name) {
            this->addGeneSegment("other", "", 0);
            bool res = this->read(filepath);
            if (is_ok) {
                *is_ok = res;
            }
        }


        GeneSegmentAlphabet(const string& name, const vector<string>& alleles, const vector<string>& sequences) {
            this->_name = name;

            this->_vec.reserve(alleles.size());

            this->addGeneSegment("other", "", 0);
            for (size_t i = 0; i < alleles.size(); i++) {
                this->addGeneSegment(alleles[i], sequences[i], i + 1);
            }
        }


        GeneSegmentAlphabet(const GeneSegmentAlphabet &other) {
            _name = other._name;
            _map.clear();
            _map.insert(other._map.begin(), other._map.end());
            _vec.clear();
            _vec.insert(_vec.begin(), other._vec.begin(), other._vec.end());
        }


        virtual ~GeneSegmentAlphabet() {}


        const string& name() const { return this->_name; }


        const GeneSegment& operator[] (segindex_t index) const { return this->_vec[index]; }
        const GeneSegment& operator[] (const string& name) const { return this->_vec[this->_map.at(name)]; }


        segindex_t size() const { return this->_vec.size(); }


        /**
        * \brief Append P(alindromic)-nucleotides to all gene segment sequences from
        * 3'end, 5'end or both.
        *
        * \param from_5_end Number of P nucleotides to append from the 3' end (i.e., sequences' start).
        * \param from_3_end Number of P nucleotides to append from the 3' end (i.e., sequences' end).
        */
        void appendPalindromicNucleotides(seq_len_t from_5_end = 4, seq_len_t from_3_end = 4) {

        }


        /**
        * \brief Save segment alphabet to the given file as a tab-separated table with 3 columns.
        */
        bool write(const string& filepath) const {
            ofstream ofs;
            ofs.open(filepath);

            if (ofs.is_open()) {
                ofs << this->_name << '\t' << "Sequences" << endl;
                for (segindex_t i = 1; i < this->_vec.size(); ++i) {
                    ofs << this->_vec[i].allele << '\t' << this->_vec[i].sequence;
                    if (i != this->_vec.size() - 1) {
                        ofs << endl;
                    }
                }
            } else {
                return false;
            }

            ofs.close();
            return true;
        }


        bool read(const string& filepath) {
            ifstream ifs;
            ifs.open(filepath);

            if (ifs.is_open()) {
                stringstream line_stream;
                string line, gene_name, gene_seq;
                bool header = false;
                segindex_t cur_index = 1;
                while (!ifs.eof()) {
                    getline(ifs, line);
                    if (line[0] != '\n') {
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
                cerr << "Gene segment alphabet [" << this->_name << "] error:" << endl << "\tinput file [" << filepath << "] not found" << endl;
                return false;
            }

            ifs.close();
            return true;
        }


    protected:

        string _name;
        unordered_map<string, segindex_t> _map;
        vector<GeneSegment> _vec;


        GeneSegmentAlphabet() {}


        void addGeneSegment(const string& name, const string& seq, segindex_t index) {
            this->_map[name] = index;
            this->_vec.push_back(GeneSegment(name, seq, index));
        }

    };


    /**
    * \class VDJRecombinationGenes
    */
    class VDJRecombinationGenes {
    public:

        VDJRecombinationGenes(const string& v_segments_name, const string& v_segments_file,
                const string& j_segments_name, const string& j_segments_file,
                bool *is_V_ok = nullptr, bool *is_J_ok = nullptr) {
            this->_V = new GeneSegmentAlphabet(v_segments_name, v_segments_file, is_V_ok);
            this->_J = new GeneSegmentAlphabet(j_segments_name, j_segments_file, is_J_ok);
            this->_D = nullptr;
        }


        VDJRecombinationGenes(const string& v_segments_name, const string& v_segments_file,
                const string& j_segments_name, const string& j_segments_file,
                const string& d_segments_name, const string& d_segments_file,
                bool *is_V_ok = nullptr, bool *is_J_ok = nullptr, bool *is_D_ok = nullptr) :
                VDJRecombinationGenes(v_segments_name, v_segments_file, v_segments_name, j_segments_file, is_V_ok, is_J_ok) {
            this->_D = new GeneSegmentAlphabet(d_segments_name, d_segments_file, is_D_ok);
        }


        VDJRecombinationGenes(const string& V_name, const vector<string>& V_alleles, const vector<string>& V_sequences,
                const string& J_name, const vector<string>& J_alleles, const vector<string>& J_sequences) {
            this->_V = new GeneSegmentAlphabet(V_name, V_alleles, V_sequences);
            this->_J = new GeneSegmentAlphabet(J_name, J_alleles, J_sequences);
            this->_D = nullptr;
        }


        VDJRecombinationGenes(const string& V_name, const vector<string>& V_alleles, const vector<string>& V_sequences,
                const string& J_name, const vector<string>& J_alleles, const vector<string>& J_sequences,
                const string& D_name, const vector<string>& D_alleles, const vector<string>& D_sequences) :
                VDJRecombinationGenes(V_name, V_alleles, V_sequences, J_name, J_alleles, J_sequences) {
            this->_D = new GeneSegmentAlphabet(D_name, D_alleles, D_sequences);
        }


        VDJRecombinationGenes(const VDJRecombinationGenes &other) {
            _V = new GeneSegmentAlphabet(*other._V);
            _J = new GeneSegmentAlphabet(*other._J);
            if (other._D) {
                _D = new GeneSegmentAlphabet(*other._D);
            } else {
                _D = nullptr;
            }
        }


        virtual ~VDJRecombinationGenes() {
            delete this->_V;
            delete this->_J;
            if (_D) { delete this->_D; }
        }


        bool is_vdj() const {
            return this->_D != nullptr;
        }


        ///@{
        const GeneSegmentAlphabet& V() const { return *(this->_V); }
        const GeneSegmentAlphabet& J() const { return *(this->_J); }
        const GeneSegmentAlphabet& D() const { return *(this->_D); }
        ///@}


        /**
        * \param segindex 0 for V(ariable), 1 for J(oining), 2 for D(iversity) gene segment.
        */
        void appendPalindromicNucleotides(GENE_SEGMENTS gene, seq_len_t from_5_end = 4, seq_len_t from_3_end = 4) {
            if (gene == VARIABLE) {
                if (this->_V) { this->_V->appendPalindromicNucleotides(from_5_end, from_3_end); }
            } else if (gene == JOINING) {
                if (this->_J) { this->_J->appendPalindromicNucleotides(from_5_end, from_3_end); }
            } else if (gene == DIVERSITY) {
                if (this->_D) { this->_D->appendPalindromicNucleotides(from_5_end, from_3_end); }
            }
        }


        bool write(const string& v_segments_file = "vsegments.txt",
                const string& j_segments_file = "jsegments.txt",
                const string& d_segments_file = "dsegments.txt") const {
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

        GeneSegmentAlphabet *_V, *_J, *_D;

        VDJRecombinationGenes() {}

    };
}

#endif