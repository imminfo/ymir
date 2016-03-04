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

#ifndef _CLONOTYPE_H
#define _CLONOTYPE_H


#include "alignment.h"


namespace ymir {


    struct Clonotype;


    /**
     * \typedef ClonotypePtr
     */
    typedef unique_ptr<Clonotype> ClonotypePtr;


    /**
    * \struct Clonotype
    */
    struct Clonotype : public VDJAlignment {

        Clonotype(const sequence_t &sequence,
                  SequenceType seq_type,
                  Recombination recomb,
                  const segments_storage_t &segments,
                  const NoGapAlignmentVector &alignments,
                  const n_D_alignments_storage_t &n_D_alignments)
                : VDJAlignment(segments, alignments, n_D_alignments),
                  _sequence(sequence),
                  _seq_type(seq_type),
                  _recomb(recomb)
        {
        }


        Clonotype(const sequence_t &sequence,
                  SequenceType seq_type,
                  Recombination recomb,
                  VDJAlignment alignment)
                : VDJAlignment(alignment),
                  _sequence(sequence),
                  _seq_type(seq_type),
                  _recomb(recomb)
        {
        }


        virtual ~Clonotype()
        {
        }


        /**
         *
         */
        ///@{
        const sequence_t& sequence() const { return _sequence; }

        const sequence_t& nuc_sequence() const {
#ifndef DNDEBUG
            check_and_throw(_seq_type != NUCLEOTIDE, "Clonotype's call to nuc_sequence() is incorrect: wrong sequence type.");
#endif
            return _sequence;
        }

        sequence_t aa_sequence() const {
#ifndef DNDEBUG
            check_and_throw(_seq_type == UNDEF_SEQ_TYPE, "Clonotype's call to aa_sequence() is incorrect: undefined sequence type.");
#endif
            if (_seq_type == NUCLEOTIDE) {
                return translate(_sequence);
            } else {
                return _sequence;
            }
        }
        ///@}


        sequence_t::const_iterator seq_iterator(seq_len_t pos) const { return _sequence.cbegin() + pos; }


        Recombination recombination() const { return _recomb; }


        SequenceType sequence_type() const { return _seq_type; }


        bool isCoding() const {
            if (_seq_type == NUCLEOTIDE) {
                return !(is_out_of_frame(_sequence) || has_end_codon(_sequence));
            } else {
                return !has_bad_aa_codons(_sequence);
            }
        }


        bool isNoncoding() const {
            if (_seq_type == NUCLEOTIDE) {
                return is_out_of_frame(_sequence) || has_end_codon(_sequence);
            } else {
                return has_bad_aa_codons(_sequence);
            }
        }


        bool isOutOfFrame() const {
            if (_seq_type == NUCLEOTIDE) {
                return is_out_of_frame(_sequence);
            } else {
                return has_oof_aa_codon(_sequence);
            }
        }

    protected:

        Recombination _recomb;

        SequenceType _seq_type;

        sequence_t _sequence; //* CDR3 or full nucleotide or amino acid sequence of a clone. */


        /**
         *
         */
        Clonotype() {}

    };

}

#endif