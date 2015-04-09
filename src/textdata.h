//
// Created by Vadim N. on 03/04/2015.
//

#ifndef YMIR_TEXTDATA_H
#define YMIR_TEXTDATA_H

#include "types.h"


namespace ymir {


    struct AbstractTDContainer {
    public:

        AbstractTDContainer(bool skip_first_column, prob_t laplace) {
            _rownames.reserve(40);
            _colnames.reserve(40);
            _laplace = laplace;
            _skip_first_column = skip_first_column;
        }

        virtual ~AbstractTDContainer() { }

        void set_laplace(prob_t laplace) { _laplace = laplace; }

        void add_row_name(const string& name) { _rownames.push_back(name); }

        void add_column_name(const string& name) { _colnames.push_back(name); }

        const vector<string>& rownames() const { return _rownames; }

        const vector<string>& colnames() const { return _colnames; }

        virtual void addRow(const vector<prob_t> vec) = 0;

        virtual bool read(const string& filepath) = 0;

    protected:

        vector<string> _rownames;
        vector<string> _colnames;
        prob_t _laplace;
        bool _skip_first_column;


        AbstractTDContainer() {
            _rownames.reserve(40);
            _colnames.reserve(40);
            _laplace = 0;
            _skip_first_column = true;
        }

    };


    struct TDVectorList : public AbstractTDContainer {
    public:

        TDVectorList(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(skip_first_column, laplace) { }

        virtual ~TDVectorList() { }

        void addRow(const vector<prob_t> vec) {

        }

        bool read(const string& filepath) {
            ifstream ifs;
            ifs.open(filepath);

            if (ifs.is_open()) {
                stringstream line_stream;
                string line, word;
                bool read_header = true;
                while (!ifs.eof()) {
                    getline(ifs, line);
                    if (line[0] != '\n') {
                        line_stream.str(line);
                        if (read_header) {
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                cout << word << endl;
                                // add row
                            }
                            read_header = false;
                        } else {
                            getline(line_stream, word, '\t'); // skip row's name
                            int i = 0;
                            while (!line_stream.eof()) {
                                getline(line_stream, word, '\t');
                                // MPFR?!?!?! I don't know
                                // add row
                                cout << word << endl;
                                ++i;
                            }
                        }
                        line_stream.clear();
                    }
                }
                return true;
            }

            return false;
        }

    protected:

        TDVectorList() : AbstractTDContainer() {}
    };


    struct TDMatrix : public AbstractTDContainer  {
    public:

        TDMatrix(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(skip_first_column, laplace)  { }

        virtual ~TDMatrix() { }

        void addRow(const vector<prob_t> vec) {

        }

        bool read(const string& filepath) {
            return false;
        }

    protected:

        TDMatrix() : AbstractTDContainer() {}

    };


    struct TDMatrixList : public AbstractTDContainer  {
    public:


        TDMatrixList(bool skip_first_column = true, prob_t laplace = 0) : AbstractTDContainer(skip_first_column, laplace)  { }

        virtual ~TDMatrixList() { }

        void addRow(const vector<prob_t> vec) {

        }

        bool read(const string& filepath) {
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
        bool ok;
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
        ok = container->read(filepath);
        if (!ok) {
            err_message  = "ERROR: problems in reading file [" + filepath + "]";
        }
        return container;
    }


    /**
    * \function read_matrix
    *
    * \brief
    *
    * \param filepath Path to the file with matrix, separated by spaces / tabs.
    *
    * \return New struct with matrix.
    */
//    NamedVectorArray read_vector_list(const string& filepath) {
//
//        ifstream ifs;
//        ifs.open(filepath);
//
//        if (ifs.is_open()) {
//            NamedVectorArray narr;
//            stringstream line_stream;
//            string line, word;
//            bool read_header = true;
//            while (!ifs.eof()) {
//                getline(ifs, line);
//                if (line[0] != '\n') {
//                    line_stream.str(line);
//                    if (read_header) {
//                        while (!line_stream.eof()) {
//                            getline(line_stream, word, '\t');
//                            narr.addColumn(word);
//                        }
//                        read_header = false;
//                    } else {
//                        getline(line_stream, word, '\t'); // skip row's name
//                        int i = 0;
//                        while (!line_stream.eof()) {
//                            getline(line_stream, word, '\t');
//                            narr.push(i, stod(word));  // MPFR?!?!?! I don't know
//                            ++i;
//                        }
//                    }
//                    line_stream.clear();
//                }
//            }
//            return narr;
//
//        } else {
//            cerr << "Matrix parsing error:" << endl << "\tinput file [" << filepath << "] not found" << endl;
//            return NamedVectorArray();
//        }
//
//        // read first line with names
//
//        // put first element from each row to column names
//
//        // remove cell at the intersection of the first row and the first column
//
//        ifs.close();
////        return NamedMatrix(colnames, rownames, Matrix<prob_t, Dynamic, Dynamic, ColMajor>(values.data()));
//    }


}

#endif //YMIR_TEXTDATA_H
