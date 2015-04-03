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
        }

        virtual ~AbstractTDContainer() { }

        void add_row_name(const string& name) { _rownames.push_back(name); }

        void add_column_name(const string& name) { _colnames.push_back(name); }

        void addRow(const vector<prob_t> vec) = 0;

        void read(const string& filepath) = 0;

    protected:

        vector<string> _rownames;
        vector<string> _colnames;

    };


    struct TDVectorList : public AbstractTDContainer {
    public:

    protected:

    };


    struct TDMatrix : public AbstractTDContainer  {
    public:

    protected:

    };


    struct TDMatrixList : public AbstractTDContainer  {
    public:

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


}

#endif //YMIR_TEXTDATA_H
