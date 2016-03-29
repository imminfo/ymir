//
// Created by Vadim N. on 24/03/2016.
//

#ifndef YMIR_CODON_TABLE_H
#define YMIR_CODON_TABLE_H


#include <unordered_map>


namespace ymir {
/*
    class CodonTable;


    class CodonTable {
    public:

//        typedef std::unordered_multimap<> aa_to_codon_storage_t;
//
//
//        typedef std::unordered_map<> codon_to_aa_storage_t;


        CodonTable()
        {
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
        }


        CodonTable(CodonTable const&) = delete;

        void operator=(CodonTable const&) = delete;


        static CodonTable& table() {
            static CodonTable table;
            return table;
        }


        std::string translate(const std::string &nuc_seq) const {
            if (nuc_seq.size() % 3 == 0 && !nuc_seq.empty()) {
                std::string res;
                // for
                return res;
            }
            return "";
        }






    private:

        aa_to_codon_storage_t _aa2codon;
        codon_to_aa_storage_t _codon2aa;

    };
*/
}

#endif //YMIR_CODON_TABLE_H
