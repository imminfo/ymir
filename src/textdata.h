//
// Created by Vadim N. on 03/04/2015.
//

#ifndef YMIR_TEXTDATA_H
#define YMIR_TEXTDATA_H

#include <vector>

#include "types.h"


namespace ymir {


    enum CONTAINER_TYPE {
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
        }


        virtual ~AbstractTDContainer() { }


        void addRowName(const string& name) { _rownames.push_back(name); }


        void addColumnName(const string& name) { _colnames.push_back(name); }


        const vector<string>& row_names() const { return _rownames; }


        const vector<string>& column_names() const { return _colnames; }


        virtual bool read(const string& filepath, string &err_message) = 0;


        const vector<prob_t>& data(size_t i) const { return _data[i]; }

        CONTAINER_TYPE type() const { return _type; }

    protected:

        vector< vector<prob_t> > _data;
        vector<string> _colnames;
        vector<string> _rownames;
        prob_t _laplace;
        bool _skip_first_column;
        CONTAINER_TYPE _type;


        AbstractTDContainer() {
            _colnames.reserve(40);
            _rownames.reserve(40);
            _data.reserve(40);
            _laplace = 0;
            _skip_first_column = true;
        }

    };


    /**
     * \brief List of vectors for deletions, insertions and gene segments.
     */
    struct TDVectorList : public AbstractTDContainer {
    public:

        TDVectorList(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(VECTOR_LIST, skip_first_column, laplace) { }


        virtual ~TDVectorList() { }


        bool read(const string& filepath, string &err_message) {
            ifstream ifs;

            ifs.open(filepath);

            if (ifs.is_open()) {
                stringstream line_stream;
                string line, word;
                vector<prob_t> word_vec;
                int line_num = 1;
                bool skip_col_num_check = true;
                while (!ifs.eof()) {
                    getline(ifs, line);
                    if (line[0] != '\n') {
                        word_vec.clear();
                        line_stream.str(line);
                        int i = !_skip_first_column;

                        if (skip_col_num_check) {
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                if (i) {
                                    this->addColumnName(word);
                                    _data.push_back(vector<prob_t>());
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
                    }

                    if (!skip_col_num_check) {
                        if (word_vec.size() == _colnames.size()) {
                            for (size_t i = 0; i < word_vec.size(); ++i) {
                                _data[i].push_back(word_vec[i]);
                            }
                        } else {
                            stringstream ss;
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

                    ++line_num;
                }
                ifs.close();
                return true;
            }

            err_message = "ERROR: can't open file [" + filepath + "]";
            return false;
        }

    protected:

        TDVectorList() : AbstractTDContainer() {}

    };


    struct TDMatrix : public AbstractTDContainer  {
    public:

        TDMatrix(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(MATRIX, skip_first_column, laplace)  { }


        virtual ~TDMatrix() { }


        bool read(const string& filepath, string &err_message) {
            ifstream ifs;

            ifs.open(filepath);

            if (ifs.is_open()) {
                _data.push_back(vector<prob_t>());
                stringstream line_stream;
                string line, word;
                vector<prob_t> word_vec;
                int line_num = 1;
                bool skip_col_num_check = true;
                while (!ifs.eof()) {
                    getline(ifs, line);
                    if (line[0] != '\n') {
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
                    }

                    if (!skip_col_num_check) {
                        if (word_vec.size() == _colnames.size()) {
                            _data[0].insert(_data[0].end(), word_vec.begin(), word_vec.end());
                        } else {
                            stringstream ss;
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

                    ++line_num;
                }
                ifs.close();
                return true;
            }

            err_message = "ERROR: can't open file [" + filepath + "]";
            return false;
        }

    protected:

        TDMatrix() : AbstractTDContainer() {}

    };


    struct TDMatrixList : public AbstractTDContainer  {
    public:


        TDMatrixList(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(MATRIX_LIST, skip_first_column, laplace)  { }


        virtual ~TDMatrixList() { }


        bool read(const string& filepath, string &err_message) {
            return false;
        }

    protected:

        TDMatrixList() : AbstractTDContainer() {}

    };


    /**
     * \function read_textdata
     *
     * \brief Read to the container data from the file with event probabilities.
     */
    AbstractTDContainer* read_textdata(const string& filepath,
                                       const string& filetype,
                                       bool skip_first_column,
                                       prob_t laplace,
                                       string& err_message) {
        AbstractTDContainer* container;
        err_message = "OK";
        if (filetype == "matrix") {
            container = new TDMatrix(skip_first_column, laplace);
        } else if (filetype == "vector.list") {
            container = new TDVectorList(skip_first_column, laplace);
        } else if (filetype == "matrix.list") {
            container = new TDMatrixList(skip_first_column, laplace);
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
