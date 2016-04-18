//
// Created by Vadim N. on 24/03/2016.
//

#ifndef YMIR_CODON_TABLE_H
#define YMIR_CODON_TABLE_H


#include <unordered_map>

#include "types.h"
#include "tools.h"


namespace ymir {

    class CodonTable;


    /**
     * \class CodonTable
     */
    class CodonTable {
    public:

        typedef std::unordered_multimap<char, sequence_t> aa_to_codon_storage_t;

        typedef std::unordered_map<int, char> codon_to_aa_storage_t;


        /**
         * \struct Codons
         */
        struct Codons {

            Codons(std::pair<aa_to_codon_storage_t::const_iterator, aa_to_codon_storage_t::const_iterator> it)
                    : _begin(it.first), _end(it.second), _current(it.first)
            {}


            bool next() {
                ++_current;
                return _current != _end;
            }


            bool end() const { return _current == _end; }


            const sequence_t& codon() const { return _current->second; }


            char aminoacid() const { return _current->first; }


            std::string print() {
                std::string res = " ";
                res[0] = this->aminoacid();
                res += ": " + this->codon() + "";
                while (this->next()) {
                    res += " ";
//                    res[res.size() - 1] = this->aminoacid();
//                    res += ":" + this->codon() + "_";
                    res += this->codon();
                }
                return res;
            }

        private:

            aa_to_codon_storage_t::const_iterator _begin, _end, _current;

            Codons() {}

        };


        /**
         *
         */
        CodonTable() {
            _aa2codon = {
                    {'A', "GCT"}, {'A', "GCC"}, {'A', "GCA"}, {'A', "GCG"},
                    {'L', "TTA"}, {'L', "TTG"}, {'L', "CTT"}, {'L', "CTC"}, {'L', "CTA"}, {'L', "CTG"},
                    {'R', "CGT"}, {'R', "CGC"}, {'R', "CGA"}, {'R', "CGG"}, {'R', "AGA"}, {'R', "AGG"},
                    {'K', "AAA"}, {'K', "AAG"},
                    {'N', "AAT"}, {'N', "AAC"},
                    {'M', "ATG"},
                    {'D', "GAT"}, {'D', "GAC"},
                    {'F', "TTT"}, {'F', "TTC"},
                    {'C', "TGT"}, {'C', "TGC"},
                    {'P', "CCT"}, {'P', "CCC"}, {'P', "CCA"}, {'P', "CCG"},
                    {'Q', "CAA"}, {'Q', "CAG"},
                    {'S', "TCT"}, {'S', "TCC"}, {'S', "TCA"}, {'S', "TCG"}, {'S', "AGT"}, {'S', "AGC"},
                    {'E', "GAA"}, {'E', "GAG"},
                    {'T', "ACT"}, {'T', "ACC"}, {'T', "ACA"}, {'T', "ACG"},
                    {'G', "GGT"}, {'G', "GGC"}, {'G', "GGA"}, {'G', "GGG"},
                    {'W', "TGG"},
                    {'H', "CAT"}, {'H', "CAC"},
                    {'Y', "TAT"}, {'Y', "TAC"},
                    {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"},
                    {'V', "GTT"}, {'V', "GTC"}, {'V', "GTA"}, {'V', "GTG"},
                    {'*', "TAA"}, {'*', "TGA"}, {'*', "TAG"}
            };

            for (aa_to_codon_storage_t::iterator it = _aa2codon.begin(); it != _aa2codon.end(); ++it) {
                sequence_t codon = it->second;
                _codon2aa[16*nuc_hash(codon[0]) + 4*nuc_hash(codon[1]) + nuc_hash(codon[2])] = it->first;
            }
        }


        CodonTable(CodonTable const&) = delete;

        void operator=(CodonTable const&) = delete;


        static CodonTable& table() {
            static CodonTable table;
            return table;
        }


        Codons codons(char aminoacid) const {
            return Codons(_aa2codon.equal_range(aminoacid));
        }


        std::string translate(const std::string &nuc_seq) const {
            if (nuc_seq.size() % 3 == 0 && !nuc_seq.empty()) {
                std::string res = "";
                for (seq_len_t i = 0; i < nuc_seq.size(); i += 3) {
                    res += _codon2aa.find(16*nuc_hash(nuc_seq[i]) + 4*nuc_hash(nuc_seq[i+1]) + nuc_hash(nuc_seq[i+2]))->second;
                }
                return res;
            }
            return "";
        }


        bool check_nucl(char aminoacid, char nucl, int pos, std::vector<bool> *vec) const {
#ifndef DNDEBUG
            assert(pos >= 0 && pos <= 2);
#endif
            auto range = _aa2codon.equal_range(aminoacid);
            vec->resize(vec->size() + 6);
            int i = 6;
            bool was_found = false, check;
            for (auto it = range.first; it != range.second; ++it) {
                check = nucl == it->second[pos];
                was_found = was_found || check;
                (*vec)[vec->size() - i] = check;
                --i;
            }
            return was_found;
        }


        codon_hash mask_nucl(char aminoacid, char nucl, int pos) const {
#ifndef DNDEBUG
            assert(pos >= 0 && pos <= 2);
#endif
            auto range = _aa2codon.equal_range(aminoacid);
            std::bitset<6> res = 0;
            int i = 0;
            for (auto it = range.first; it != range.second; ++it) {
                res.set(5 - i, nucl == it->second[pos]);
                ++i;
            }
            return res.to_ulong();
        }


        std::string print() {
            std::string res = "";
            for (aa_to_codon_storage_t::iterator it = _aa2codon.begin(); it != _aa2codon.end(); ++it) {
                res += this->codons(it->first).print();
                res += '\n';
            }
            return res;
        }

        const aa_to_codon_storage_t aminoacids() const { return _aa2codon; }

    private:

        aa_to_codon_storage_t _aa2codon;
        codon_to_aa_storage_t _codon2aa;

    };


    /**
     * \function translate
     *
     * \brief Translate the given nucleotide sequence.
     */
    sequence_t translate(const sequence_t& nuc_seq) {
        return CodonTable::table().translate(nuc_seq);
    }

}

#endif //YMIR_CODON_TABLE_H
