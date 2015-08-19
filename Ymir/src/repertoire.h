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

#ifndef _REPERTOIRE_H_
#define _REPERTOIRE_H_


#include <memory>

#include "aligner.h"
#include "clonotype.h"
#include "maag.h"
#include "tools.h"


namespace ymir {

    #define DEFAULT_REPERTOIRE_RESERVE_SIZE 300000


    typedef std::vector<MAAG> MAAGRepertoire;


    class ClonesetView {

    public:

        ClonesetView() : _source(new std::vector<Clonotype>()) {
            _shifts.resize(0);
        }


        ClonesetView(std::shared_ptr<std::vector<Clonotype>> pvec, const std::vector<size_t>& shifts) :
                _source(pvec), _shifts(shifts) {}


        ClonesetView(const ClonesetView& other) : _source(other._source), _shifts(other._shifts) {}


        virtual ~ClonesetView() { }


        const Clonotype& operator[] (size_t index) const {
            return _source->at(_shifts[index]);
        }
        ClonesetView subvec(const std::vector<size_t> &indices) const {
            std::vector<size_t> shifts;
            shifts.reserve(indices.size());
            for (size_t i = 0; i < indices.size(); ++i) {
                shifts.push_back(_shifts[indices[i]]);
            }
            return ClonesetView(_source, shifts);
        };


        size_t size() const { return _shifts.size(); }


        ClonesetView coding() const {
            std::vector<size_t> inds;
            inds.reserve(this->size() / 3);
            for (size_t i = 0; i < this->size(); ++i) {
                if (!is_out_of_frame((*this)[i].sequence()) && !has_end_codon((*this)[i].sequence())) {
                    inds.push_back(i);
                }
            }
            return this->subvec(inds);
        }


        ClonesetView noncoding(bool out_of_frames_only = false) const {
            std::vector<size_t> inds;
            inds.reserve(this->size() / 3);
            for (size_t i = 0; i < this->size(); ++i) {
                if (is_out_of_frame((*this)[i].sequence()) || (!out_of_frames_only && has_end_codon((*this)[i].sequence()))) {
                    inds.push_back(i);
                }
            }
            return this->subvec(inds);
        }


        ClonesetView head(size_t size) const {
            std::vector<size_t> shifts;
            shifts.reserve(size);
            for (size_t i = 0; i < size; ++i) {
                shifts.push_back(i);
            }
            return this->subvec(shifts);
        };


        ClonesetView slice(size_t start, size_t end) const {
            if (end > start) {
                std::vector<size_t> shifts;
                shifts.reserve(end - start + 1);
                for (size_t i = start; i < end; ++i) {
                    shifts.push_back(i);
                }
                return this->subvec(shifts);
            } else {
                cerr << "Error in ClonesetView: the end index is lower than the start index." << endl;
                return ClonesetView();
            }
        };


    protected:

        std::shared_ptr<std::vector<Clonotype>> _source;
        std::vector<size_t> _shifts;

    };


    class Cloneset : public ClonesetView {

    public:


        Cloneset() : ClonesetView() {
            _shifts.resize(0);
        }


        // swap constructor
        Cloneset(std::vector<Clonotype>& vec) {
            this->swap(vec);
        }


        virtual ~Cloneset() { }


        void swap(std::vector<Clonotype>& vec) {
            this->_source->swap(vec);
            this->_shifts.resize(this->_source->size());
            for (size_t i = 0; i < this->_source->size(); ++i) {
                this->_shifts[i] = i;
            }
        }


        Clonotype& operator[] (size_t index) {
            return _source->at(_shifts[index]);
        }


    protected:

    };

}

#endif