//
// Created by Vadim N. on 03/04/2015.
//


// Disclaimer: I doesn't like IO in C++.


#ifndef YMIR_TEXTDATA_H
#define YMIR_TEXTDATA_H

//#include <vector>

#include "types.h"


namespace ymir {


    enum CONTAINER_TYPE {
        VECTOR,
        VECTOR_LIST,
        MATRIX,
        MATRIX_LIST
    };


    struct AbstractTDContainer {
    public:

        AbstractTDContainer(CONTAINER_TYPE ctype, bool skip_first_column, prob_t laplace) {
            _colnames.reserve(40);
            _rownames.reserve(40);
            _data.reserve(40);
            _laplace = laplace;
            _skip_first_column = skip_first_column;
            _type = ctype;
            _file_exists = false;
        }


        virtual ~AbstractTDContainer() { }


        virtual void addDataValue(prob_t value, size_t i = 0) { _data[i].push_back(value + _laplace); }


        virtual void addDataVector(const std::vector<prob_t> &vec) {
            _data.push_back(std::vector<prob_t>());
            for (auto i = 0; i < vec.size(); ++i) {
                _data[_data.size() - 1].push_back(vec[i] + _laplace);
            }
        }
        
        
        virtual void addDataVector(std::vector<prob_t>::const_iterator start, std::vector<prob_t>::const_iterator end) {
            _data.push_back(std::vector<prob_t>());
            _data[_data.size() - 1].insert(_data[_data.size() - 1].cend(), start, end);
        }


        prob_t laplace() const { return _laplace; }


        seg_index_t n_rows() const { return _rownames.size(); }


        seg_index_t n_columns() const { return _colnames.size(); }


        void addRowName(const std::string& name) { _rownames.push_back(name); }


        void addColumnName(const std::string& name) { _colnames.push_back(name); }


        void addMetadata(seq_len_t data) { _metadata.push_back(data); }


        const std::vector<std::string>& row_names() const { return _rownames; }


        const std::vector<std::string>& column_names() const { return _colnames; }


        seq_len_t metadata(size_t i) const { return _metadata[i]; }


        virtual bool read(const std::string& filepath, std::string &err_message) = 0;


        virtual bool write(const std::string& filepath) = 0;


        const std::vector<prob_t>& data(size_t i) const { return _data[i]; }


        CONTAINER_TYPE type() const { return _type; }


        bool file_exists() const { return _file_exists; }

    protected:

        std::vector< std::vector<prob_t> > _data;
        std::vector<std::string> _colnames;
        std::vector<std::string> _rownames;
        prob_t _laplace;
        bool _skip_first_column, _file_exists;
        CONTAINER_TYPE _type;
        std::vector<seq_len_t> _metadata;


        AbstractTDContainer() {
            _colnames.reserve(40);
            _rownames.reserve(40);
            _data.reserve(40);
            _laplace = 0;
            _skip_first_column = true;
        }

    };


    /**
     * \brief Vector of gene segments.
     */
    struct TDVector : public AbstractTDContainer {
    public:

        TDVector(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(VECTOR, skip_first_column, laplace) { }


        virtual ~TDVector() { }


        bool read(const std::string& filepath, std::string &err_message) {
            ifstream ifs;

            ifs.open(filepath);

            if (ifs.is_open()) {
                std::stringstream line_stream;
                std::string line, word;
                std::vector<prob_t> word_vec;
                int line_num = 1;
                bool skip_col_num_check = true, empty_line = false;

                while (!ifs.eof()) {
                    getline(ifs, line);
                    if (line[0] != '\n' && line[0] != 0) {
                        word_vec.clear();
                        line_stream.str(line);
                        int i = !_skip_first_column;

                        if (skip_col_num_check) {
                            getline(line_stream, word, '\t');
                            this->addColumnName("Probability");
                            _data.push_back(std::vector<prob_t>());
                            _data[0].reserve(60);
                        } else {
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                // MPFR?!?!?! I don't know
                                // add row
                                if (i) { word_vec.push_back(stod(word)); }
                                else { this->addRowName(word); }
                                ++i;
                            }

                        }
                        line_stream.clear();
                    } else {
                        empty_line = true;
                    }

                    if (!empty_line) {
                        if (!skip_col_num_check) {
                            if (word_vec.size() == _colnames.size()) {
                                for (size_t i = 0; i < word_vec.size(); ++i) {
                                    _data[i].push_back(word_vec[i]);
                                }
                            } else {
                                std::stringstream ss;
                                ss << "ERROR: number of elements doesn't match the number of columns in the line " <<
                                (int) line_num <<
                                " (expected: " <<
                                (int) _colnames.size() <<
                                ", got: " <<
                                (int) word_vec.size() << ")";
                                err_message = ss.str();
                                return false;
                            }
                        } else {
                            skip_col_num_check = false;
                        }
                    } else {
                        empty_line = false;
                    }

                    ++line_num;
                }
                ifs.close();
                _file_exists = true;
                return true;
            }

            _file_exists = false;
            err_message = "ERROR: can't open file [" + filepath + "]";
            return false;
        }


        bool write(const std::string& filepath) {
            ofstream ofs;

            ofs.open(filepath);

            if (ofs.is_open()) {

                ofs << _colnames[0] << '\t' << _colnames[1] << std::endl;

                for (auto i = 0; i < _data[0].size(); ++i) {
                    ofs << _rownames[i] << '\t' << _data[0][i] << std::endl;
                }

                ofs.close();
                _file_exists = true;
                return true;
            }

            _file_exists = false;
            return false;
        }

    protected:

        TDVector() : AbstractTDContainer() {}

    };



    /**
     * \brief List of std::vectors for deletions and insertions.
     */
    struct TDVectorList : public AbstractTDContainer {
    public:

        TDVectorList(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(VECTOR_LIST, skip_first_column, laplace) { }


        virtual ~TDVectorList() { }


        bool read(const std::string& filepath, std::string &err_message) {
            ifstream ifs;

            ifs.open(filepath);

            if (ifs.is_open()) {
                std::stringstream line_stream;
                std::string line, word;
                std::vector<prob_t> word_vec;
                int line_num = 1;
                bool skip_col_num_check = true, empty_line = false;
                while (!ifs.eof()) {
                    getline(ifs, line);
                    if (line[0] != '\n' && line[0] != 0) {
                        word_vec.clear();
                        line_stream.str(line);
                        int i = !_skip_first_column;

                        if (skip_col_num_check) {
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                if (i) {
                                    this->addColumnName(word);
                                    _data.push_back(std::vector<prob_t>());
                                    _data[_data.size() - 1].reserve(40);
                                }
                                ++i;
                            }
                        } else {
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                // MPFR?!?!?! I don't know
                                // add row
                                if (i) { word_vec.push_back(stod(word)); }
                                else { this->addRowName(word); }
                                ++i;
                            }

                        }
                        line_stream.clear();
                    } else {
                        empty_line = true;
                    }

                    if (!empty_line) {
                        if (!skip_col_num_check) {
                            if (word_vec.size() == _colnames.size()) {
                                for (size_t i = 0; i < word_vec.size(); ++i) {
                                    _data[i].push_back(word_vec[i]);
                                }
                            } else {
                                std::stringstream ss;
                                ss << "ERROR: number of elements doesn't match the number of columns in the line " <<
                                (int) line_num <<
                                " (expected: " <<
                                (int) _colnames.size() <<
                                ", got: " <<
                                (int) word_vec.size() << ")";
                                err_message = ss.str();
                                return false;
                            }
                        } else {
                            skip_col_num_check = false;
                        }
                    } else {
                        empty_line = false;
                    }

                    ++line_num;
                }
                ifs.close();
                _file_exists = true;
                return true;
            }

            _file_exists = false;
            err_message = "ERROR: can't open file [" + filepath + "]";
            return false;
        }


        bool write(const std::string& filepath) {
            ofstream ofs;

            ofs.open(filepath);

            if (ofs.is_open()) {
                seq_len_t max_len = _data[0].size();
                for (auto i = 1; i < _data.size(); ++i) { max_len = std::max(max_len, (seq_len_t) _data[i].size()); }

                for (auto i = 0; i < _colnames.size(); ++i) {
                    ofs << _colnames[i] << "\t";
                }
                ofs << std::endl;

                for (auto pos_i = 0; pos_i < max_len; ++pos_i) {
                    if (pos_i < _rownames.size()) {
                        ofs << _rownames[pos_i] << '\t';
                    } else {
                        ofs << "0" << '\t';
                    }

                    for (auto vec_i = 0; vec_i < _data.size(); ++vec_i) {
                        if (pos_i < _data[vec_i].size()) {
                            ofs << _data[vec_i][pos_i];
                        } else {
                            ofs << 0;
                        }
                        if (vec_i < _data.size() - 1) { ofs << '\t'; }
                    }
                    ofs << std::endl;
                }

                ofs.close();
                _file_exists = true;
                return true;
            }

            _file_exists = false;
            return false;
        }

    protected:

        TDVectorList() : AbstractTDContainer() {}

    };


    struct TDMatrix : public AbstractTDContainer  {
    public:

        TDMatrix(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(MATRIX, skip_first_column, laplace)  { }


        virtual ~TDMatrix() { }


        virtual void addDataVector(const std::vector<prob_t> &vec) {
            if (!_data.size()) {
                _data.push_back(std::vector<prob_t>());
            }
            _data[0].insert(_data[0].end(), vec.begin(), vec.end());
        }


        bool read(const std::string& filepath, std::string &err_message) {
            ifstream ifs;

            ifs.open(filepath);

            if (ifs.is_open()) {
                _data.push_back(std::vector<prob_t>());
                std::stringstream line_stream;
                std::string line, word;
                std::vector<prob_t> word_vec;
                int line_num = 1;
                bool skip_col_num_check = true, empty_line = false;
                while (!ifs.eof()) {
                    getline(ifs, line);
                    if (line[0] != '\n' && line[0] != 0) {
                        word_vec.clear();
                        line_stream.str(line);
                        int i = !_skip_first_column;

                        if (skip_col_num_check) {
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                if (i) { this->addColumnName(word); }
                                ++i;
                            }
                        } else {
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                // MPFR?!?!?! I don't know
                                // add row
                                if (i) { word_vec.push_back(stod(word)); }
                                else { this->addRowName(word); }
                                ++i;
                            }

                        }
                        line_stream.clear();
                    } else {
                        empty_line = true;
                    }

                    if (!empty_line) {
                        if (!skip_col_num_check) {
                            if (word_vec.size() == _colnames.size()) {
                                _data[0].insert(_data[0].end(), word_vec.begin(), word_vec.end());
                            } else {
                                std::stringstream ss;
                                ss << "ERROR: number of elements doesn't match the number of columns in the line " <<
                                (int) line_num <<
                                " (expected: " <<
                                (int) _colnames.size() <<
                                ", got: " <<
                                (int) word_vec.size() << ")";
                                err_message = ss.str();
                                return false;
                            }
                        } else {
                            skip_col_num_check = false;
                        }
                    } else {
                        empty_line = false;
                    }

                    ++line_num;
                }
                ifs.close();
                _file_exists = true;
                return true;
            }

            _file_exists = false;
            err_message = "ERROR: can't open file [" + filepath + "]";
            return false;
        }


        bool write(const std::string& filepath) {
            ofstream ofs;

            ofs.open(filepath);

            if (ofs.is_open()) {
                for (auto i = 0; i < _colnames.size(); ++i) {
                    ofs << _colnames[i] << "\t";
                }
                ofs << std::endl;

                for (size_t row_i = 0; row_i < _data[0].size() / _metadata[0]; ++row_i) {
                    ofs << _rownames[row_i] << '\t';
                    for (auto col_i = 0; col_i < _metadata[0]; ++col_i) {
                        ofs << _data[0][row_i * _metadata[0] + col_i];
                        if (col_i < _metadata[0] - 1) { ofs << '\t'; }
                    }
                    ofs << std::endl;
                }

                ofs.close();
                _file_exists = true;
                return true;
            }

            _file_exists = false;
            return false;
        }

    protected:

        TDMatrix() : AbstractTDContainer() {}

    };


    /**
     * row names - segments
     * col names - no col names
     * columns - number of columns in each matrices
     */
    struct TDMatrixList : public AbstractTDContainer  {
    public:


        TDMatrixList(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(MATRIX_LIST, skip_first_column, laplace)  { }


        virtual ~TDMatrixList() { }


        bool read(const std::string& filepath, std::string &err_message) {
            ifstream ifs;

            ifs.open(filepath);

            if (ifs.is_open()) {
                std::stringstream line_stream;
                std::string line, word;
                std::vector<prob_t> word_vec;
                int line_num = 1, matrix_num = 0;
                bool skip_col_num_check = true, in_matrix = true, empty_line = false;
                _data.push_back(std::vector<prob_t>());
                while (!ifs.eof()) {
                    getline(ifs, line);
                    if (line[0] != '\n' && line[0] != 0) {
                        word_vec.clear();
                        line_stream.str(line);
                        int i = !_skip_first_column;

                        if (skip_col_num_check) {
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                if (i) { this->addColumnName(word); }
                                else { this->addRowName(word); }
                                ++i;
                            }
                            _metadata.push_back(i - 1);
                        } else {
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                // MPFR?!?!?! I don't know
                                // add row
                                if (i) { word_vec.push_back(stod(word)); }
                                ++i;
                            }

                        }
                        line_stream.clear();
                    } else {
                        matrix_num++;
                        _data.push_back(std::vector<prob_t>());
                        empty_line = true;
                        skip_col_num_check = true;
                    }

                    if (!empty_line) {
                        if (!skip_col_num_check) {
                            if (word_vec.size() == _metadata[matrix_num]) {
                                _data[matrix_num].insert(_data[matrix_num].end(), word_vec.begin(), word_vec.end());
                            } else {
                                std::stringstream ss;
                                ss << "ERROR: number of elements doesn't match the number of columns in the line " <<
                                (int) line_num <<
                                " (expected: " <<
                                (int) _colnames.size() <<
                                ", got: " <<
                                (int) word_vec.size() << ")";
                                err_message = ss.str();
                                return false;
                            }
                        } else {
                            skip_col_num_check = false;
                        }
                    } else {
                        empty_line = false;
                    }

                    ++line_num;
                }
                ifs.close();
                _file_exists = true;
                return true;
            }

            _file_exists = false;
            err_message = "ERROR: can't open file [" + filepath + "]";
            return false;
        }


        bool write(const std::string& filepath) {
            ;
        }

    protected:

        TDMatrixList() : AbstractTDContainer() {}

    };


    /**
     * \function read_textdata
     *
     * \brief Read to the container data from the file with event probabilities.
     */
    AbstractTDContainer* read_textdata(const std::string& filepath,
                                       const std::string& filetype,
                                       bool skip_first_column,
                                       prob_t laplace,
                                       std::string& err_message) {
        AbstractTDContainer* container;
        err_message = "OK";
        if (filetype == "matrix") {
            container = new TDMatrix(skip_first_column, laplace);
        } else if (filetype == "vector.list") {
            container = new TDVectorList(skip_first_column, laplace);
        } else if (filetype == "matrix.list") {
            container = new TDMatrixList(skip_first_column, laplace);
        } else if (filetype == "vector") {
            container = new TDVector(skip_first_column, laplace);
        } else {
            if (filetype == "") {
                err_message = "ERROR: no file type for [" + filepath + "]";
            } else {
                err_message = "ERROR: unrecognised file type for [" + filepath + "]";
            }
            return nullptr;
        }

        container->read(filepath, err_message);
        return container;
    }

}

#endif //YMIR_TEXTDATA_H
