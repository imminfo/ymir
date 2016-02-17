//
// Created by Vadim N. on 14/02/2016.
//

#ifndef YMIR_ALIGNMENT_MATRIX_H
#define YMIR_ALIGNMENT_MATRIX_H


#include <algorithm>
#include <string>
#include <vector>

#include "nogap_alignment_vector.h"
#include "gapped_alignment_vector.h"
#include "types.h"


namespace ymir {

    /**
     * \struct AlignmentMatrixBase
     */
    template <typename AlignmentType>
    class AlignmentMatrixBase {

    public:

        typedef std::vector<bool> bit_storage_t;


        typedef std::vector<alignment_score_t > score_storage_t;


        static const seq_len_t default_nrows = 80;


        static const seq_len_t default_ncols = 20;


        AlignmentMatrixBase()
                : _gene(0),
                  _nrow(default_nrows),
                  _ncol(default_ncols),
                  _starts(_nrow * _ncol, false),
                  _matrix(_nrow * _ncol, 0)
        {
        }


        AlignmentMatrixBase(seg_index_t gene,
                            const sequence_t &pattern,
                            const sequence_t &text)
                : _gene(gene),
                  _nrow(pattern.size() + 1),
                  _ncol(text.size() + 1),
                  _starts(_nrow * _ncol, false),
                  _matrix(_nrow * _ncol, 0)
        {
        }


        void reinit(seg_index_t gene,
                    const sequence_t &pattern,
                    const sequence_t &text)
        {
            _gene = gene;
            _nrow = pattern.size() + 1;
            _ncol = text.size() + 1;
            _starts.resize(_nrow * _ncol);
            std::fill(_starts.begin(), _starts.end(), false);
            _matrix.resize(_nrow * _ncol);
            std::fill(_matrix.begin(), _matrix.end(), 0);
        }


        void setStart(seq_len_t row, seq_len_t col) { _starts[index(row, col)] = true; }


        alignment_score_t getBestAlignment(AlignmentType *vec, const sequence_t &pattern, const sequence_t &text) const;


        /**
         *
         */
        ///@{
        alignment_score_t score(seq_len_t row, seq_len_t col) const { return _matrix[index(row, col)]; }

        alignment_score_t& score(seq_len_t row, seq_len_t col) { return _matrix[index(row, col)]; }
        ///@}

    private:

        seq_len_t _nrow, _ncol;
        seg_index_t _gene;
        bit_storage_t _starts;
        score_storage_t _matrix;


        size_t index(seq_len_t row, seq_len_t col) const { return row * _ncol + col; }

    };


    typedef AlignmentMatrixBase<NoGapAlignmentVector> SWNGAlignmentMatrix;


    typedef AlignmentMatrixBase<GappedAlignmentVector> SWAlignmentMatrix;


    template <>
    alignment_score_t SWNGAlignmentMatrix::getBestAlignment(NoGapAlignmentVector *vec,
                                                            const sequence_t &pattern,
                                                            const sequence_t &text) const
    {
        AlignmentVectorBase::events_storage_t bitvec;

        // Find maximum score.
        seq_len_t max_i = 0, max_j = 0;
        alignment_score_t max_score = -1;
        for (seq_len_t i = 0; i < _nrow; ++i) {
            for (seq_len_t j = 0; j < _ncol; ++j) {
                if (score(i, j) > max_score) {
                    max_i = i;
                    max_j = j;
                    max_score = score(i, j);
                }
            }
        }

        // Traceback to the start, storing alignment events.
        seq_len_t cur_i = max_i, cur_j = max_j;
        while (!_starts[index(cur_i, cur_j)]) {
            bitvec.push_back(pattern[cur_i] == text[cur_j]);
            --cur_i;
            --cur_j;
        }

        vec->addAlignment(_gene, cur_i, cur_j, bitvec);

        return score(max_i, max_j);
    }


    template <>
    alignment_score_t SWAlignmentMatrix::getBestAlignment(GappedAlignmentVector *vec,
                                                          const sequence_t &pattern,
                                                          const sequence_t &text) const
    {
        AlignmentVectorBase::events_storage_t bitvec;

        // Find maximum score.
        seq_len_t max_i = 0, max_j = 0;
        alignment_score_t max_score = -1;
        for (seq_len_t i = 0; i < _nrow; ++i) {
            for (seq_len_t j = 0; j < _ncol; ++j) {
                if (score(i, j) > max_score) {
                    max_i = i;
                    max_j = j;
                    max_score = score(i, j);
                }
            }
        }

        // Traceback to the start, storing alignment events.
        seq_len_t cur_i = max_i, cur_j = max_j, max_index = 0;
        std::array<alignment_score_t, 3> score_arr;
        while (!_starts[index(cur_i, cur_j)]) {
            score_arr[0] = score(cur_i - 1, cur_j - 1);
            score_arr[1] = cur_j > 0 ? score(cur_i, cur_j - 1) : -1;
            score_arr[2] = cur_i > 0 ? score(cur_i - 1, cur_j) : -1;
            max_index = std::distance(score_arr.begin(), std::max_element(score_arr.begin(), score_arr.end()));
            switch (max_index) {
                case 0:
//                    std::cout  << (int) cur_i << ":" << (int) cur_j << std::endl;
                    pattern[cur_i - 1] == text[cur_j - 1] ? add_match(&bitvec) : add_mismatch(&bitvec);
                    --cur_i;
                    --cur_j;
                    break;

                case 1:
                    add_ins(&bitvec);
                    --cur_j;
                    break;

                case 2:
                    add_del(&bitvec);
                    --cur_i;
                    break;

                default:
                    break;
            }
        }

        add_match(&bitvec);
        // reverse
        AlignmentVectorBase::events_storage_t bitvec2;

        for (int i = bitvec.size() / 2 - 1; i >= 0; --i) {
            bitvec2.push_back(bitvec[i*2]);
            bitvec2.push_back(bitvec[i*2 + 1]);
        }
        vec->addAlignment(_gene, cur_i, cur_j, bitvec2);

        return score(max_i, max_j);
    }

}

#endif //YMIR_ALIGNMENT_MATRIX_H
