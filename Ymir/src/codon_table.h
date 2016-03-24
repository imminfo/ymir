//
// Created by Vadim N. on 24/03/2016.
//

#ifndef YMIR_CODON_TABLE_H
#define YMIR_CODON_TABLE_H


#include <unordered_map>


namespace ymir {

    class CodonTable;


    class CodonTable {
    public:

        typedef std::unordered_multimap<> aa_to_codon_storage_t;


        typedef std::unordered_map<> codon_to_aa_storage_t;


        CodonTable()
        {
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

}

#endif //YMIR_CODON_TABLE_H
