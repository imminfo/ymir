//
// Created by Vadim N. on 07/08/2015.
//

#ifndef YMIR_WRITER_H
#define YMIR_WRITER_H


#include "repertoire.h"


namespace ymir {

    class RepertoireWriter {
    public:


        RepertoireWriter() { }


        virtual ~RepertoireWriter() { }


        bool write(const std::string &filepath, const ClonesetView &cloneset) const {

        }

    protected:


    };
}


#endif //YMIR_WRITER_H
