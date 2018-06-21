//
// Created by Vadim N. on 09/04/2016.
//

#ifndef YMIR_MAAG_BASE_H
#define YMIR_MAAG_BASE_H


#include "multimatrixchain3.h"

namespace ymir {

#define VJ_CHAIN_SIZE 4
#define VDJ_CHAIN_SIZE 7


    /**
     * \class MAAGBase
     */
    class MAAGBase : protected ProbMMC {

        friend class MAAGBuilder;
        friend class MAAGForwardBackwardAlgorithm;

    public:


        /**
         * \brief Default constructor.
         */
        MAAGBase()
            : _recomb(UNDEF_RECOMB),
              _sequence(nullptr),
              _seq_poses(nullptr),
              _n_poses(0),
              _seq_type(UNDEF_SEQ_TYPE)
        {
            _values.reserve(1);
        }


        MAAGBase(SequenceType seq_type)
            : _recomb(UNDEF_RECOMB),
              _sequence(nullptr),
              _seq_poses(nullptr),
              _n_poses(0),
              _seq_type(seq_type)
        {
            _values.reserve(1);
        }


        /**
         *
         */
        MAAGBase(const MAAGBase &other)
            : ProbMMC(other),
              _recomb(other._recomb),
              _n_poses(other._n_poses),
              _seq_type(other._seq_type)
        {
            if (other._seq_poses) {
                _seq_poses.reset(new seq_len_t[_n_poses]);
                std::copy(other._seq_poses.get(), other._seq_poses.get() + _n_poses, _seq_poses.get());
            } else {
                _seq_poses = nullptr;
            }

            if (other._sequence) { _sequence.reset(new std::string(*other._sequence)); }
            else { _sequence = nullptr; }
        }


        MAAGBase(MAAGBase &&other) {
            std::swap(_recomb, other._recomb);

            _chain.swap(other._chain);
            _values.swap(other._values);

            std::swap(_n_poses, other._n_poses);
            _seq_poses.swap(other._seq_poses);

            _sequence.swap(other._sequence);
            std::swap(_seq_type, other._seq_type);
        }


        /**
         *
         */
        virtual ~MAAGBase() { }


        MAAGBase& operator= (const MAAGBase &other) {
            _recomb = other._recomb;
            _chain = other._chain;
            _values = other._values;

            _n_poses = other._n_poses;

            if (other._seq_poses) {
                _seq_poses.reset(new seq_len_t[_n_poses]);
                std::copy(other._seq_poses.get(), other._seq_poses.get() + _n_poses, _seq_poses.get());
            } else {
                _seq_poses = nullptr;
            }

            if (other._sequence) {
                _sequence.reset(new std::string(*other._sequence));
            }
            else { _sequence.reset(); }

            _seq_type = other._seq_type;

            return *this;
        }


        MAAGBase& operator=(MAAGBase &&other) {
            std::swap(_recomb, other._recomb);

            _chain.swap(other._chain);
            _values.swap(other._values);

            std::swap(_n_poses, other._n_poses);
            _seq_poses.swap(other._seq_poses);

            _sequence.swap(other._sequence);
            std::swap(_seq_type, other._seq_type);

            return *this;
        }


        /**
         * \brief Get the number of aligned gene segments.
         *
         * \return Number of the aligned specific gene segments.
         */
        ///@{
        event_ind_t nVar() const { return (_recomb == VJ_RECOMB) ? this->nodeRows(VJ_VAR_JOI_GEN_I) : this->nodeSize(VDJ_VAR_GEN_I); }

        event_ind_t nJoi() const { return (_recomb == VJ_RECOMB) ? this->nodeColumns(VJ_VAR_JOI_GEN_I) : this->nodeSize(VDJ_JOI_DEL_I); }

        event_ind_t nDiv() const { return (_recomb == VJ_RECOMB) ? 0 : _chain[VDJ_DIV_DEL_I].size(); }
        ///@}


        /**
         *
         */
        seq_len_t position(seq_len_t i) const { return _seq_poses[i]; }


        /**
         * \brief Save the graph to the harddrive and load a previously saved graph.
         *
         * \param filepath Path to the file with a saved graph for loading and new file for saving.
         * \param stream An output / input stream, from that read the graph or save the graph.
         */
        ///@{
        virtual bool save(const std::string &filepath) const {
            std::ofstream ofs(filepath);
            return this->save(ofs);
        }

        virtual bool save(std::ostream &stream) const {
            return false;
        }

        virtual bool load(const std::string &filepath) {
            std::ifstream ifs(filepath);
            return this->load(ifs);
        }

        virtual bool load(std::istream &stream) {
            return false;
        }
        ///@}


        /**
         *
         */
        ///@{
        Recombination recombination() const { return _recomb; }

        bool is_vj() const { return _recomb == VJ_RECOMB; }

        bool is_vdj() const { return _recomb == VDJ_RECOMB; }

//        bool is_vd2j() const { return _recomb == VD2J_RECOMB; }
        ///@}


        const std::string& sequence() const {
#ifndef DNDEBUG
            if (!_sequence) { throw(std::runtime_error("Access to a MAAG sequence when it's a nullptr!")); }
#endif
            return *_sequence;
        }


        seq_len_t n_poses() const { return _n_poses; }


        SequenceType sequence_type() const { return _seq_type; }

        dim_t rows(node_ind_t node_i) const { return this->nodeRows(node_i); }

        dim_t cols(node_ind_t node_i) const { return this->nodeColumns(node_i); }


    protected:

        unique_ptr<seq_len_t[]> _seq_poses;  /** Vector of the initial clonotype sequence's positions for each vertex. 0-based */
        unique_ptr<sequence_t> _sequence;  /** Nucleotide or amino acid CDR3 sequence. */

        SequenceType _seq_type;
        Recombination _recomb;

        seq_len_t _n_poses;

    };
}

#endif //YMIR_MAAG_BASE_H
