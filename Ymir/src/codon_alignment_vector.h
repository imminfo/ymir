//
// Created by Vadim N. on 11/04/2016.
//

#ifndef YMIR_CODON_ALIGNMENT_VECTOR_H
#define YMIR_CODON_ALIGNMENT_VECTOR_H


#include "alignment_vector_base.h"


namespace ymir {


    typedef uint8_t codon_hash;


    codon_hash compute_codon_hash(const AlignmentVectorBase::events_storage_t &bits, seq_len_t start) {
        return codon_hash(  (bits[start]     << 5)
                          + (bits[start + 1] << 4)
                          + (bits[start + 2] << 3)
                          + (bits[start + 3] << 2)
                          + (bits[start + 4] << 1)
                          + (bits[start + 5]));
    }


    struct CodonAlignmentVector : public AlignmentVectorBase {


        CodonAlignmentVector() : AlignmentVectorBase()
        {
        }


        void addAlignment(seg_index_t id, seq_len_t p_start, seq_len_t t_start, const events_storage_t &vec) {
#ifndef DNDEBUG
            assert(vec.size() % 6 == 0);
#endif

            _data.push_back(p_start);
            _data.push_back(t_start);
            _data.push_back(vec.size() / 6);
            _data.push_back(id);
            _starts.push_back(_events.size());
            _events.insert(_events.end(), vec.begin(), vec.end());
        }


        /**
         * \param i 0-based index of pattern.
         * \param j 1-based position.
         */
        codon_hash getCodon(seq_len_t i, seq_len_t j) const {
#ifndef DNDEBUG
            check_and_throw(_events.size() == 0, "Codon alignment vector: no codons stored in the vector.");

            if (_starts[i] + 6 * (j - 1) >= _events.size()) {
                std::cout << (int) (_starts[i] + 6 * (j - 1)) << std::endl;
                std::cout << (int) _starts.size() << std::endl;
                std::cout << (int) i << std::endl;
                std::cout << (int) _starts[i] << std::endl;
                std::cout << (int) j << std::endl;
                std::cout << (int) _events.size() << std::endl;
                throw(std::runtime_error("Codon alignment vector: codon index is out of bounds."));
            }
#endif
            return compute_codon_hash(_events, _starts[i] + 6 * (j - 1));
        }

    };

}

#endif //YMIR_CODON_ALIGNMENT_VECTOR_H
