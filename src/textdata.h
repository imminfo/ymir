//
// Created by Vadim N. on 03/04/2015.
//

#ifndef YMIR_TEXTDATA_H
#define YMIR_TEXTDATA_H

#include "types.h"


namespace ymir {


    struct AbstractTDContainer {
    public:

        AbstractTDContainer() {
            _rownames.reserve(40);
            _colnames.reserve(40);
            _laplace = 0;
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

    };


    struct TDVectorList : public AbstractTDContainer {
    public:

        TDVectorList() : AbstractTDContainer() { }

        virtual ~TDVectorList() { }

        void addRow(const vector<prob_t> vec) {

        }

        bool read(const string& filepath) {
            return false;
        }

    protected:

    };


    struct TDMatrix : public AbstractTDContainer  {
    public:

        TDMatrix() : AbstractTDContainer()  { }

        virtual ~TDMatrix() { }

        void addRow(const vector<prob_t> vec) {

        }

        bool read(const string& filepath) {
            return false;
        }

    protected:

    };


    struct TDMatrixList : public AbstractTDContainer  {
    public:


        TDMatrixList() : AbstractTDContainer()  { }

        virtual ~TDMatrixList() { }

        void addRow(const vector<prob_t> vec) {

        }

        bool read(const string& filepath) {
            return false;
        }

    protected:

    };


    /**
     * \struct NamedVectoArray
     */
    struct NamedVectorArray {

    public:

        struct Vector {

        public:

            vector<prob_t> vec;
            string name;

            Vector(const string& name_) {
                name = name_;
                vec.clear();
                vec.reserve(10);
            }
        };

        NamedVectorArray() {
            cols.clear();
            cols.reserve(10);
        }

        void addColumn(const string& name) {
            cols.push_back(name);
            map[name] = cols.size() - 1;
        }

        Vector& operator[](const string& name) { return cols[map[name]]; }

        Vector& operator[](int index) { return cols[index]; }

        void push(int index, prob_t value) { cols[index].vec.push_back(value); }

        string getName(int index) const { return cols[index].name; }


        vector<Vector> cols;
        unordered_map<string, int> map;
    };


    /**
    * \function read_matrix
    *
    * \brief
    *
    * \param filepath Path to the file with matrix, separated by spaces / tabs.
    *
    * \return New struct with matrix.
    */
    NamedVectorArray read_vector_list(const string& filepath) {

        ifstream ifs;
        ifs.open(filepath);

        if (ifs.is_open()) {
            NamedVectorArray narr;
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
                            narr.addColumn(word);
                        }
                        read_header = false;
                    } else {
                        getline(line_stream, word, '\t'); // skip row's name
                        int i = 0;
                        while (!line_stream.eof()) {
                            getline(line_stream, word, '\t');
                            narr.push(i, stod(word));  // MPFR?!?!?! I don't know
                            ++i;
                        }
                    }
                    line_stream.clear();
                }
            }
            return narr;

        } else {
            cerr << "Matrix parsing error:" << endl << "\tinput file [" << filepath << "] not found" << endl;
            return NamedVectorArray();
        }

        // read first line with names

        // put first element from each row to column names

        // remove cell at the intersection of the first row and the first column

        ifs.close();
//        return NamedMatrix(colnames, rownames, Matrix<prob_t, Dynamic, Dynamic, ColMajor>(values.data()));
    }


}

#endif //YMIR_TEXTDATA_H
