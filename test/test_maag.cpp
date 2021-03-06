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
//#define NCURSES_TERM_H_incl 1


#include "testutils.h"


using namespace std;
using namespace ymir;


std::string TEST_DATA_FOLDER;


YMIR_TEST_START(test_mmc)
    /*
    0-1:
        .5
    0-2:
        2

    1-1:
        1 2 4
    1-2:
        2 3 0

    2:
        1 1 1
        2 2 0
        1 0 .5

    3-1:
        2
        4
        2
    3-2:
        3
        4
        5

    4-1:
        .2
    4-2:
        .4

    prod 1-1-1 = 4.4
    prod 2-1-2 = 52.8
     */

    ProbMMC mat;
    YMIR_ASSERT(mat.addNode(2, 1, 1) == 0)
    mat(0, 0, 0, 0) = .5;
    mat(0, 1, 0, 0) = 2;
    YMIR_ASSERT(mat.nodeRows(0) == 1)
    YMIR_ASSERT(mat.nodeColumns(0) == 1)

    YMIR_ASSERT(mat.addNode(2, 1, 3) == 1)
    mat(1, 0, 0, 0) = 1;
    mat(1, 0, 0, 1) = 2;
    mat(1, 0, 0, 2) = 4;
    mat(1, 1, 0, 0) = 2;
    mat(1, 1, 0, 1) = 3;
    mat(1, 1, 0, 2) = 0;
    YMIR_ASSERT(mat.nodeRows(1) == 1)
    YMIR_ASSERT(mat.nodeColumns(1) == 3)

    YMIR_ASSERT(mat.addNode(1, 3, 3) == 2)
    mat(2, 0, 0, 0) = 1;
    mat(2, 0, 0, 1) = 1;
    mat(2, 0, 0, 2) = 1;
    mat(2, 0, 1, 0) = 2;
    mat(2, 0, 1, 1) = 2;
    mat(2, 0, 1, 2) = 0;
    mat(2, 0, 2, 0) = 1;
    mat(2, 0, 2, 1) = 0;
    mat(2, 0, 2, 2) = .5;
    YMIR_ASSERT(mat.nodeRows(2) == 3)
    YMIR_ASSERT(mat.nodeColumns(2) == 3)

    YMIR_ASSERT(mat.addNode(2, 3, 1) == 3)
    mat(3, 0, 0, 0) = 2;
    mat(3, 0, 1, 0) = 4;
    mat(3, 0, 2, 0) = 2;
    mat(3, 1, 0, 0) = 3;
    mat(3, 1, 1, 0) = 4;
    mat(3, 1, 2, 0) = 5;

    YMIR_ASSERT(mat.addNode(2, 1, 1) == 4)
    mat(4, 0, 0, 0) = .2;
    mat(4, 1, 0, 0) = .4;

    YMIR_ASSERT((mat.matrix(0, 0)
            * mat.matrix(1, 0)
            * mat.matrix(2, 0)
            * mat.matrix(3, 0)
            * mat.matrix(4, 0))(0, 0) - 4.4 < 1e-14)

    YMIR_ASSERT((mat.matrix(0, 1)
            * mat.matrix(1, 1)
            * mat.matrix(2, 0)
            * mat.matrix(3, 1)
            * mat.matrix(4, 1))(0, 0) - 52.8 < 1e-14)

YMIR_TEST_END


YMIR_TEST_START(test_maag_vj)

    ModelParameterVector mvec = make_test_events_vj();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("CCCG");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCCGAG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
    seqvec2.push_back("CCGTTT");
    seqvec2.push_back("ATTT");
    seqvec2.push_back("AGGTTT");

    VDJRecombinationGenes genes("VA", alvec1, seqvec1, "JA", alvec2, seqvec2);

    MAAGBuilder maag_builder(mvec, genes);

    ClonotypeNucBuilder cl_builder;
    // CCCGACGGTTT
    // .....CCGTTT
    // .......ATTT
    // .....AGGTTT
    // CCCG.AC.GGTTT
    // poses:
    // vs: 0-1-2-3-4-5
    // js: 7-8-9-10-11-12
    cl_builder.setSequence("CCCGACGGTTT")
            .addVarAlignment(1, 1, 1, 4)
            .addVarAlignment(3, 1, 1, 5)
            .addJoiAlignment(1, 3, 8, 4)
            .addJoiAlignment(2, 2, 9, 3)
            .addJoiAlignment(3, 2, 7, 5);
    cl_builder.setRecombination(VJ_RECOMB);
    ClonotypeNuc clonotype = cl_builder.buildClonotype();

//    cout << "here" << endl;
    MAAGnuc maag = maag_builder.build(clonotype, SAVE_METADATA, NO_ERRORS);

    YMIR_ASSERT2(maag.nVar(), 2)
    YMIR_ASSERT2(maag.nJoi(), 3)
    YMIR_ASSERT2(maag.nDiv(), 0)

    YMIR_ASSERT2(maag.position(0), 0)
    YMIR_ASSERT2(maag.position(3), 3)
    YMIR_ASSERT2(maag.position(5), 5)
    YMIR_ASSERT2(maag.position(6), 7)
    YMIR_ASSERT2(maag.position(8), 9)
    YMIR_ASSERT2(maag.position(10), 11)
    YMIR_ASSERT2(maag.position(11), 12)

//    cout << "VJ 0:" << maag.rows(0) << ":" << maag.cols(0) << endl;
//    cout << "Vdel 1:" << maag.rows(1) << ":" << maag.cols(1) << endl;
//    cout << "VJins 2:" << maag.rows(2) << ":" << maag.cols(2) << endl;
//    cout << "Jdel 3:" << maag.rows(3) << ":" << maag.cols(3) << endl;

    YMIR_ASSERT2(maag.event_index(0, 0, 0, 0), mvec.event_index(VJ_VAR_JOI_GEN, 0, 0, 0))
    YMIR_ASSERT2(maag.event_index(0, 0, 0, 1), mvec.event_index(VJ_VAR_JOI_GEN, 0, 0, 1))
    YMIR_ASSERT2(maag.event_index(0, 0, 1, 0), mvec.event_index(VJ_VAR_JOI_GEN, 0, 2, 0))
    YMIR_ASSERT2(maag.event_index(0, 0, 1, 2), mvec.event_index(VJ_VAR_JOI_GEN, 0, 2, 2))

    YMIR_ASSERT2(maag.event_index(1, 0, 0, 0), mvec.event_index(VJ_VAR_DEL, 0, 4));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 1), mvec.event_index(VJ_VAR_DEL, 0, 3));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 2), mvec.event_index(VJ_VAR_DEL, 0, 2));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 3), mvec.event_index(VJ_VAR_DEL, 0, 1));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 4), mvec.event_index(VJ_VAR_DEL, 0, 0));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 5), 0)

    YMIR_ASSERT2(maag.event_index(1, 1, 0, 0), mvec.event_index(VJ_VAR_DEL, 2, 6))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 1), mvec.event_index(VJ_VAR_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 2), mvec.event_index(VJ_VAR_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 3), mvec.event_index(VJ_VAR_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 4), mvec.event_index(VJ_VAR_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 5), mvec.event_index(VJ_VAR_DEL, 2, 1))

    YMIR_ASSERT2(maag.event_index(2, 0, 0, 0), mvec.event_index(VJ_VAR_JOI_INS_LEN, 0, maag.position(6) - maag.position(0) - 1))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 1), 0)

    YMIR_ASSERT2(maag.event_index(3, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 1, 0), mvec.event_index(VJ_JOI_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(3, 0, 2, 0), mvec.event_index(VJ_JOI_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(3, 0, 3, 0), mvec.event_index(VJ_JOI_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(3, 0, 4, 0), mvec.event_index(VJ_JOI_DEL, 0, 5))
    YMIR_ASSERT2(maag.event_index(3, 0, 5, 0), mvec.event_index(VJ_JOI_DEL, 0, 6))

    YMIR_ASSERT2(maag.event_index(3, 1, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 1, 1, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 1, 2, 0), mvec.event_index(VJ_JOI_DEL, 1, 1))
    YMIR_ASSERT2(maag.event_index(3, 1, 3, 0), mvec.event_index(VJ_JOI_DEL, 1, 2))
    YMIR_ASSERT2(maag.event_index(3, 1, 4, 0), mvec.event_index(VJ_JOI_DEL, 1, 3))
    YMIR_ASSERT2(maag.event_index(3, 1, 5, 0), mvec.event_index(VJ_JOI_DEL, 1, 4))

    YMIR_ASSERT2(maag.event_index(3, 2, 0, 0), mvec.event_index(VJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_index(3, 2, 1, 0), mvec.event_index(VJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(3, 2, 2, 0), mvec.event_index(VJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(3, 2, 3, 0), mvec.event_index(VJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(3, 2, 4, 0), mvec.event_index(VJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(3, 2, 5, 0), mvec.event_index(VJ_JOI_DEL, 2, 6))

YMIR_TEST_END


YMIR_TEST_START(test_maag_vj_err)

    ModelParameterVector mvec = make_test_events_vj();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("CCCG");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCCGAG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
    seqvec2.push_back("CCGTTT");
    seqvec2.push_back("ATTT");
    seqvec2.push_back("AGGTTT");

    VDJRecombinationGenes genes("VA", alvec1, seqvec1, "JA", alvec2, seqvec2);

    MAAGBuilder maag_builder(mvec, genes);

    NoGapAlignmentVector vec;
    NoGapAlignmentVector::events_storage_t bits;
    ClonotypeNucBuilder cl_builder;
    // CCCGACGGTTT
    // CCCG
    // CCCGAG
    // .....CCGTTT
    // .......ATTT
    // .....AGGTTT
    // CCCG.AC.GGTTT
    // poses:
    // vs: 0-1-2-3-4-5
    // js: 7-8-9-10-11-12
    cl_builder.setSequence("CCCGACGGTTT").setRecombination(VJ_RECOMB);

    bits = {0, 0, 0, 0};
    vec.addAlignment(1, 1, 1, bits);

    bits = {0, 0, 0, 0, 0, 1};
    vec.addAlignment(3, 1, 1, bits);

    cl_builder.addVarAlignment(vec);
    vec.clear();

    bits = {0, 1, 0, 0, 0, 0};
    vec.addAlignment(1, 1, 6, bits);

    bits = {1, 0, 0, 0};
    vec.addAlignment(2, 1, 8, bits);

    bits = {1, 0, 0, 0, 0, 0};
    vec.addAlignment(3, 1, 6, bits);

    cl_builder.addJoiAlignment(vec);
    vec.clear();

    ClonotypeNuc clonotype = cl_builder.buildClonotype();

//    cout << "here" << endl;
    MAAGnuc maag = maag_builder.build(clonotype, SAVE_METADATA, COMPUTE_ERRORS);

    YMIR_ASSERT2(maag.nVar(), 2)
    YMIR_ASSERT2(maag.nJoi(), 3)
    YMIR_ASSERT2(maag.nDiv(), 0)

    YMIR_ASSERT(clonotype.is_good())
    YMIR_ASSERT(maag.has_events())
    YMIR_ASSERT(maag.has_errors())

    YMIR_ASSERT2(maag.position(0), 0)
    YMIR_ASSERT2(maag.position(3), 3)
    YMIR_ASSERT2(maag.position(5), 5)
    YMIR_ASSERT2(maag.position(6), 6)
    YMIR_ASSERT2(maag.position(8), 7)
    YMIR_ASSERT2(maag.position(10), 11)
    YMIR_ASSERT2(maag.position(11), 12)

    YMIR_ASSERT2(maag.event_index(0, 0, 0, 0), mvec.event_index(VJ_VAR_JOI_GEN, 0, 0, 0))
    YMIR_ASSERT2(maag.event_index(0, 0, 0, 1), mvec.event_index(VJ_VAR_JOI_GEN, 0, 0, 1))
    YMIR_ASSERT2(maag.event_index(0, 0, 1, 0), mvec.event_index(VJ_VAR_JOI_GEN, 0, 2, 0))
    YMIR_ASSERT2(maag.event_index(0, 0, 1, 2), mvec.event_index(VJ_VAR_JOI_GEN, 0, 2, 2))
    YMIR_ASSERT2(maag.errors(0, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.errors(0, 0, 0, 1), 0)
    YMIR_ASSERT2(maag.errors(0, 0, 0, 2), 0)
    YMIR_ASSERT2(maag.errors(0, 0, 0, 3), 0)
    YMIR_ASSERT2(maag.errors(0, 0, 0, 4), 0)
    YMIR_ASSERT2(maag.errors(0, 0, 0, 5), 0)
    YMIR_ASSERT2(maag.errors(0, 0, 0, 6), 0)

    YMIR_ASSERT2(maag.event_index(1, 0, 0, 0), mvec.event_index(VJ_VAR_DEL, 0, 4));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 1), mvec.event_index(VJ_VAR_DEL, 0, 3));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 2), mvec.event_index(VJ_VAR_DEL, 0, 2));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 3), mvec.event_index(VJ_VAR_DEL, 0, 1));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 4), mvec.event_index(VJ_VAR_DEL, 0, 0));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 5), 0)
    YMIR_ASSERT2(maag.errors(0, 1, 0, 0), 0)
    YMIR_ASSERT2(maag.errors(0, 1, 0, 1), 0)
    YMIR_ASSERT2(maag.errors(0, 1, 0, 2), 0)
    YMIR_ASSERT2(maag.errors(0, 1, 0, 3), 0)
    YMIR_ASSERT2(maag.errors(0, 1, 0, 4), 0)
    YMIR_ASSERT2(maag.errors(0, 1, 0, 5), 0)
    YMIR_ASSERT2(maag.errors(0, 1, 0, 6), 1)

    YMIR_ASSERT2(maag.event_index(1, 1, 0, 0), mvec.event_index(VJ_VAR_DEL, 2, 6))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 1), mvec.event_index(VJ_VAR_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 2), mvec.event_index(VJ_VAR_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 3), mvec.event_index(VJ_VAR_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 4), mvec.event_index(VJ_VAR_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 5), mvec.event_index(VJ_VAR_DEL, 2, 1))

    YMIR_ASSERT2(maag.event_index(2, 0, 0, 0), mvec.event_index(VJ_VAR_JOI_INS_LEN, 0, maag.position(6) - maag.position(0) - 1))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 1), 0)

    YMIR_ASSERT2(maag.event_index(3, 0, 0, 0), mvec.event_index(VJ_JOI_DEL, 0, 0))
    YMIR_ASSERT2(maag.event_index(3, 0, 1, 0), mvec.event_index(VJ_JOI_DEL, 0, 1))
    YMIR_ASSERT2(maag.event_index(3, 0, 2, 0), mvec.event_index(VJ_JOI_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(3, 0, 3, 0), mvec.event_index(VJ_JOI_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(3, 0, 4, 0), mvec.event_index(VJ_JOI_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(3, 0, 5, 0), mvec.event_index(VJ_JOI_DEL, 0, 5))
    YMIR_ASSERT2(maag.event_index(3, 0, 6, 0), mvec.event_index(VJ_JOI_DEL, 0, 6))
    YMIR_ASSERT2(maag.errors(1, 0, 0, 0), 1)
    YMIR_ASSERT2(maag.errors(1, 0, 1, 0), 1)
    YMIR_ASSERT2(maag.errors(1, 0, 2, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 0, 3, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 0, 4, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 0, 5, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 0, 6, 0), 0)

    YMIR_ASSERT2(maag.event_index(3, 1, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 1, 1, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 1, 2, 0), mvec.event_index(VJ_JOI_DEL, 1, 0))
    YMIR_ASSERT2(maag.event_index(3, 1, 3, 0), mvec.event_index(VJ_JOI_DEL, 1, 1))
    YMIR_ASSERT2(maag.event_index(3, 1, 4, 0), mvec.event_index(VJ_JOI_DEL, 1, 2))
    YMIR_ASSERT2(maag.event_index(3, 1, 5, 0), mvec.event_index(VJ_JOI_DEL, 1, 3))
    YMIR_ASSERT2(maag.event_index(3, 1, 6, 0), mvec.event_index(VJ_JOI_DEL, 1, 4))
    YMIR_ASSERT2(maag.errors(1, 1, 0, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 1, 1, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 1, 2, 0), 1)
    YMIR_ASSERT2(maag.errors(1, 1, 3, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 1, 4, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 1, 5, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 1, 6, 0), 0)

    YMIR_ASSERT2(maag.event_index(3, 2, 0, 0), mvec.event_index(VJ_JOI_DEL, 2, 0))
    YMIR_ASSERT2(maag.event_index(3, 2, 1, 0), mvec.event_index(VJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_index(3, 2, 2, 0), mvec.event_index(VJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(3, 2, 3, 0), mvec.event_index(VJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(3, 2, 4, 0), mvec.event_index(VJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(3, 2, 5, 0), mvec.event_index(VJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(3, 2, 6, 0), mvec.event_index(VJ_JOI_DEL, 2, 6))
    YMIR_ASSERT2(maag.errors(1, 2, 0, 0), 1)
    YMIR_ASSERT2(maag.errors(1, 2, 1, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 2, 2, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 2, 3, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 2, 4, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 2, 5, 0), 0)
    YMIR_ASSERT2(maag.errors(1, 2, 6, 0), 0)

YMIR_TEST_END


YMIR_TEST_START(test_maag_vdj)

    ModelParameterVector mvec = make_test_events_vdj();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("CCCG");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCCGGG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
    seqvec2.push_back("CCGTTT");
    seqvec2.push_back("ATTT");
    seqvec2.push_back("AGGTTT");

    vector<string> alvec3;
    vector<string> seqvec3;
    alvec3.push_back("Dseg1");
    alvec3.push_back("Dseg2");
    alvec3.push_back("Dseg3");
    seqvec3.push_back("GTTT");
    seqvec3.push_back("ACCGT");
    seqvec3.push_back("ACCACC");

    VDJRecombinationGenes genes("VB", alvec1, seqvec1, "JB", alvec2, seqvec2, "DB", alvec3, seqvec3);

    MAAGBuilder maag_builder(mvec, genes);

    ClonotypeNucBuilder cl_builder;

    /*
     D1:
     CCCGACCGGTTT
             GTTT

     D2:
     CCCGACCGGTTT
         ACCG.T
    A.CCG.T

     D3:
       CCCGACCGGTTT
           ACC.ACC
       ACC.ACC
    */
    cl_builder.setSequence("CCCGACCGGTTT")
            .addVarAlignment(1, 1, 1, 4)
            .addVarAlignment(3, 1, 1, 5)

            .addJoiAlignment(1, 3, 9, 4)
            .addJoiAlignment(2, 2, 10, 3)
            .addJoiAlignment(3, 2, 8, 5)

            .addDivAlignment(2, 1, 5, 4)
            .addDivAlignment(2, 2, 2, 3)
            .addDivAlignment(3, 1, 5, 3)
            .addDivAlignment(3, 4, 5, 3)
            .addDivAlignment(1, 1, 9, 4);
    cl_builder.setRecombination(VDJ_RECOMB);
    ClonotypeNuc clonotype = cl_builder.buildClonotype();

    MAAGnuc maag = maag_builder.build(clonotype, SAVE_METADATA, NO_ERRORS);

    cout << "V :" << maag.rows(0) << ":" << maag.cols(0) << endl;
    cout << "Vdel :" << maag.rows(1) << ":" << maag.cols(1) << endl;
    cout << "VDins :" << maag.rows(2) << ":" << maag.cols(2) << endl;
    cout << "Ddel :" << maag.rows(3) << ":" << maag.cols(3) << endl;
    cout << "DJins :" << maag.rows(4) << ":" << maag.cols(4) << endl;
    cout << "Jdel :" << maag.rows(5) << ":" << maag.cols(5) << endl;
    cout << "JD :" << maag.rows(6) << ":" << maag.cols(6) << endl;

    YMIR_ASSERT2(maag.nVar(), 2)
    YMIR_ASSERT2(maag.nJoi(), 3)
    YMIR_ASSERT2(maag.nDiv(), 3)

    YMIR_ASSERT2(maag.event_index(0, 0, 0, 0), mvec.event_index(VDJ_VAR_GEN, 0, 0))
    YMIR_ASSERT2(maag.event_index(0, 1, 0, 0), mvec.event_index(VDJ_VAR_GEN, 0, 2))

    YMIR_ASSERT2(maag.event_index(1, 0, 0, 0), mvec.event_index(VDJ_VAR_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 1), mvec.event_index(VDJ_VAR_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 2), mvec.event_index(VDJ_VAR_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 3), mvec.event_index(VDJ_VAR_DEL, 0, 1))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 4), mvec.event_index(VDJ_VAR_DEL, 0, 0))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 5), 0)

    YMIR_ASSERT2(maag.event_index(1, 1, 0, 0), mvec.event_index(VDJ_VAR_DEL, 2, 6))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 1), mvec.event_index(VDJ_VAR_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 2), mvec.event_index(VDJ_VAR_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 3), mvec.event_index(VDJ_VAR_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 4), mvec.event_index(VDJ_VAR_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 5), mvec.event_index(VDJ_VAR_DEL, 2, 1))

//    // seq: 1 2 4 6 7 8 9
//    YMIR_ASSERT2(maag.event_index(2, 0, 0, 0), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 0))
//    YMIR_ASSERT2(maag.event_index(2, 0, 0, 1), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 1))
//    YMIR_ASSERT2(maag.event_index(2, 0, 0, 2), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 3))
//    YMIR_ASSERT2(maag.event_index(2, 0, 0, 3), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 5))
//    YMIR_ASSERT2(maag.event_index(2, 0, 2, 2), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 1))
//
//    rows: 2 5 5 6 9 10
//    // rows 1 2 3 4 5 | cols 3 4 5 6 7
//    // TODO: check D deletions indices
//    YMIR_ASSERT2(maag.event_index(3, 0, 0, 0), 0)
//    YMIR_ASSERT2(maag.event_index(3, 0, 0, 1), 0)
//    YMIR_ASSERT2(maag.event_index(3, 0, 1, 0), 0)
//    YMIR_ASSERT2(maag.event_index(3, 0, 1, 1), mvec.event_index(VDJ_DIV_DEL, 1, 2, 2))
//    YMIR_ASSERT2(maag.event_index(3, 0, 2, 2), 0)
//    YMIR_ASSERT2(maag.event_index(3, 0, 3, 3), mvec.event_index(VDJ_DIV_DEL, 1, 3, 1))
//    YMIR_ASSERT2(maag.event_index(3, 0, 3, 4), mvec.event_index(VDJ_DIV_DEL, 1, 3, 0))
//    YMIR_ASSERT2(maag.event_index(3, 0, 4, 4), mvec.event_index(VDJ_DIV_DEL, 1, 4, 0))
//
//    YMIR_ASSERT2(maag.event_index(3, 2, 5, 5), mvec.event_index(VDJ_DIV_DEL, 0, 0, 3))
//    YMIR_ASSERT2(maag.event_index(3, 2, 6, 5), 0)
//    YMIR_ASSERT2(maag.event_index(3, 2, 5, 6), mvec.event_index(VDJ_DIV_DEL, 0, 1, 0))
//
//    // row: 3 4 6 7 8 10 11
//    // col: 7 8 9 10 11 [12]
//    YMIR_ASSERT2(maag.event_index(4, 0, 0, 0), mvec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 3))
//    YMIR_ASSERT2(maag.event_index(4, 0, 1, 1), mvec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 3))
//    YMIR_ASSERT2(maag.event_index(4, 0, 3, 0), 0)

    YMIR_ASSERT2(maag.event_index(5, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(5, 0, 1, 0), mvec.event_index(VDJ_JOI_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(5, 0, 2, 0), mvec.event_index(VDJ_JOI_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(5, 0, 3, 0), mvec.event_index(VDJ_JOI_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(5, 0, 4, 0), mvec.event_index(VDJ_JOI_DEL, 0, 5))
    YMIR_ASSERT2(maag.event_index(5, 0, 5, 0), mvec.event_index(VDJ_JOI_DEL, 0, 6))

    YMIR_ASSERT2(maag.event_index(5, 1, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(5, 1, 1, 0), 0)
    YMIR_ASSERT2(maag.event_index(5, 1, 2, 0), mvec.event_index(VDJ_JOI_DEL, 1, 1))
    YMIR_ASSERT2(maag.event_index(5, 1, 3, 0), mvec.event_index(VDJ_JOI_DEL, 1, 2))
    YMIR_ASSERT2(maag.event_index(5, 1, 4, 0), mvec.event_index(VDJ_JOI_DEL, 1, 3))
    YMIR_ASSERT2(maag.event_index(5, 1, 5, 0), mvec.event_index(VDJ_JOI_DEL, 1, 4))

    YMIR_ASSERT2(maag.event_index(5, 2, 0, 0), mvec.event_index(VDJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_index(5, 2, 1, 0), mvec.event_index(VDJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(5, 2, 2, 0), mvec.event_index(VDJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(5, 2, 3, 0), mvec.event_index(VDJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(5, 2, 4, 0), mvec.event_index(VDJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(5, 2, 5, 0), mvec.event_index(VDJ_JOI_DEL, 2, 6))

    YMIR_ASSERT2(maag.event_index(6, 0, 0, 0), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 0, 1))
    YMIR_ASSERT2(maag.event_index(6, 0, 0, 1), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 0, 2))
    YMIR_ASSERT2(maag.event_index(6, 0, 0, 2), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 0, 0))
    YMIR_ASSERT2(maag.event_index(6, 0, 2, 0), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 2, 1))
    YMIR_ASSERT2(maag.event_index(6, 0, 2, 1), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 2, 2))
    YMIR_ASSERT2(maag.event_index(6, 0, 2, 2), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 2, 0))

//    cout << "save" << endl;
//    cout << maag_builder.build(clonotype, SAVE_METADATA).fullProbability() << endl;
//    cout << "no-save" << endl;
//    cout << maag_builder.build(clonotype, NO_METADATA).fullProbability() << endl;

    // i don't want to compute by hand all this crazy matrices!
    // i've already tested chain products in previous tests!
    // ):<
    // also it's (A, B) not (A - B < eps) because this results are pretty precise on this toy example
    YMIR_ASSERT2(maag.fullProbability(0, 0, 0), maag_builder.build(clonotype, NO_METADATA, NO_ERRORS).fullProbability(0, 0, 0))  // error is here
    YMIR_ASSERT2(maag.fullProbability(1, 1, 1), maag_builder.build(clonotype, NO_METADATA, NO_ERRORS).fullProbability(1, 1, 1))
    YMIR_ASSERT2(maag.fullProbability(0, 2, 2), maag_builder.build(clonotype, NO_METADATA, NO_ERRORS).fullProbability(0, 2, 2))  // error is here

YMIR_TEST_END


YMIR_TEST_START(test_maag_vdj_err)

    ModelParameterVector mvec = make_test_events_vdj();

    YMIR_ASSERT(mvec.error_prob() != 0)

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("CCCG");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCCGGG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
    seqvec2.push_back("CCGTTT");
    seqvec2.push_back("ATTT");
    seqvec2.push_back("AGGTTT");

    vector<string> alvec3;
    vector<string> seqvec3;
    alvec3.push_back("Dseg1");
    alvec3.push_back("Dseg2");
    alvec3.push_back("Dseg3");
    seqvec3.push_back("GTTT");
    seqvec3.push_back("ACCGGT");
    seqvec3.push_back("CCCGGAC");

    VDJRecombinationGenes genes("VB", alvec1, seqvec1, "JB", alvec2, seqvec2, "DB", alvec3, seqvec3);

    MAAGBuilder maag_builder(mvec, genes);

    NoGapAlignmentVector vec, vec2;
    NoGapAlignmentVector::events_storage_t bits;
    ClonotypeNucBuilder cl_builder;

    cl_builder.setSequence("CCCGACGGTTT").setRecombination(VDJ_RECOMB);

    bits = {0, 0, 0, 0};
    vec.addAlignment(1, 1, 1, bits);

    bits = {0, 0, 0, 0, 0, 1};
    vec.addAlignment(3, 1, 1, bits);

    cl_builder.addVarAlignment(vec);
    vec.clear();

    bits = {0, 1, 0, 0, 0, 0};
    vec.addAlignment(1, 1, 6, bits);

    bits = {1, 0, 0, 0};
    vec.addAlignment(2, 1, 8, bits);
    bits = {1, 0, 0, 0, 0, 0};
    vec.addAlignment(3, 1, 6, bits);

    cl_builder.addJoiAlignment(vec);
    vec.clear();

    /*
     D1:
            CCGTTT
             |||||
            AGGTTT
       CCCGACGGTTT
             .GTTT
             .GTT.T
             G.TTT
     D2:
       CCCGACGGTTT
       ACCGGT
          ACCGGT

     D3:
       CCCGACGGTTT
      CCCGGAC
       CCCGGAC
          CCCGGAC
    */
    bits = {1, 0, 0, 0, 1, 1};
    vec.addAlignment(2, 1, 1, bits);
    bits = {1, 1, 0, 0, 0, 0};
    vec2.addAlignment(2, 1, 4, bits);
    vec.extend(vec2); vec2.clear();
    cl_builder.addDivAlignment(vec);
    vec.clear();

    bits = {0, 0, 1, 0, 0, 0};
    vec.addAlignment(3, 2, 1, bits);
    bits = {0, 0, 0, 0, 1, 1, 1};
    vec2.addAlignment(3, 1, 1, bits);
    vec.extend(vec2); vec2.clear();
    bits = {1, 1, 0, 0, 0, 1, 1};
    vec2.addAlignment(3, 1, 4, bits);
    vec.extend(vec2); vec2.clear();
    cl_builder.addDivAlignment(vec);
    vec.clear();

    bits = {0, 0, 0, 0};
    vec.addAlignment(1, 1, 8, bits);
    cl_builder.addDivAlignment(vec);

    ClonotypeNuc clonotype = cl_builder.buildClonotype();

    MAAGnuc maag = maag_builder.build(clonotype, SAVE_METADATA, COMPUTE_ERRORS);

    YMIR_ASSERT2(maag.nVar(), 2)
    YMIR_ASSERT2(maag.nJoi(), 3)
    YMIR_ASSERT2(maag.nDiv(), 3)

    YMIR_ASSERT2(maag.event_index(0, 0, 0, 0), mvec.event_index(VDJ_VAR_GEN, 0, 0))
    YMIR_ASSERT2(maag.event_index(0, 1, 0, 0), mvec.event_index(VDJ_VAR_GEN, 0, 2))

    YMIR_ASSERT2(maag.event_index(1, 0, 0, 0), mvec.event_index(VDJ_VAR_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 1), mvec.event_index(VDJ_VAR_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 2), mvec.event_index(VDJ_VAR_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 3), mvec.event_index(VDJ_VAR_DEL, 0, 1))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 4), mvec.event_index(VDJ_VAR_DEL, 0, 0))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 5), 0)

    YMIR_ASSERT2(maag.event_index(1, 1, 0, 0), mvec.event_index(VDJ_VAR_DEL, 2, 6))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 1), mvec.event_index(VDJ_VAR_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 2), mvec.event_index(VDJ_VAR_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 3), mvec.event_index(VDJ_VAR_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 4), mvec.event_index(VDJ_VAR_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 5), mvec.event_index(VDJ_VAR_DEL, 2, 1))

    // seq: 1 2 4 6 7 8 9
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 0), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 0))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 1), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 1))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 2), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 3))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 3), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 5))
    YMIR_ASSERT2(maag.event_index(2, 0, 2, 2), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 1))

    // rows 1 2 3 4 5 6 7 8 | cols
    // TODO: check D deletions indices
    YMIR_ASSERT2(maag.event_index(3, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 0, 1), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 1, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 1, 1), mvec.event_index(VDJ_DIV_DEL, 1, 2, 2))
    YMIR_ASSERT2(maag.event_index(3, 0, 2, 2), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 3, 3), mvec.event_index(VDJ_DIV_DEL, 1, 3, 1))
    YMIR_ASSERT2(maag.event_index(3, 0, 3, 4), mvec.event_index(VDJ_DIV_DEL, 1, 3, 0))
    YMIR_ASSERT2(maag.event_index(3, 0, 4, 4), mvec.event_index(VDJ_DIV_DEL, 1, 4, 0))
    // bits = {1, 0, 0, 0, 1, 1};
    YMIR_ASSERT2(maag.errors(1, 0, 0, 0), 1)
    YMIR_ASSERT2(maag.errors(1, 0, 0, 1), 1)
    YMIR_ASSERT2(maag.errors(1, 0, 0, 2), 2)
    YMIR_ASSERT2(maag.errors(1, 0, 0, 3), 3)
    YMIR_ASSERT2(maag.errors(1, 0, 1, 3), 2)
    YMIR_ASSERT2(maag.errors(1, 0, 2, 4), 1)
    YMIR_ASSERT2(maag.errors(1, 0, 3, 5), 2)

    YMIR_ASSERT2(maag.errors(1, 1, 0, 0), 0)

    YMIR_ASSERT2(maag.event_index(3, 2, 5, 5), mvec.event_index(VDJ_DIV_DEL, 0, 0, 3))
    YMIR_ASSERT2(maag.event_index(3, 2, 6, 5), 0)
    YMIR_ASSERT2(maag.event_index(3, 2, 5, 6), mvec.event_index(VDJ_DIV_DEL, 0, 1, 0))
    YMIR_ASSERT2(maag.errors(1, 2, 0, 0), 1)

    // row: 3 4 6 7 8 10 11
    // col: 7 8 9 10 11 [12]
    YMIR_ASSERT2(maag.event_index(4, 0, 0, 0), mvec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 3))
    YMIR_ASSERT2(maag.event_index(4, 0, 1, 1), mvec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 3))
    YMIR_ASSERT2(maag.event_index(4, 0, 3, 0), 0)

    YMIR_ASSERT2(maag.event_index(5, 0, 0, 0), mvec.event_index(VDJ_JOI_DEL, 0, 0))
    YMIR_ASSERT2(maag.event_index(5, 0, 1, 0), mvec.event_index(VDJ_JOI_DEL, 0, 1))
    YMIR_ASSERT2(maag.event_index(5, 0, 2, 0), mvec.event_index(VDJ_JOI_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(5, 0, 3, 0), mvec.event_index(VDJ_JOI_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(5, 0, 4, 0), mvec.event_index(VDJ_JOI_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(5, 0, 5, 0), mvec.event_index(VDJ_JOI_DEL, 0, 5))
    YMIR_ASSERT2(maag.event_index(5, 0, 6, 0), mvec.event_index(VDJ_JOI_DEL, 0, 6))
    YMIR_ASSERT2(maag.errors(2, 0, 0, 0), 1)
    YMIR_ASSERT2(maag.errors(2, 0, 1, 0), 1)
    YMIR_ASSERT2(maag.errors(2, 0, 2, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 0, 3, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 0, 4, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 0, 5, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 0, 6, 0), 0)

    YMIR_ASSERT2(maag.event_index(5, 1, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(5, 1, 1, 0), 0)
    YMIR_ASSERT2(maag.event_index(5, 1, 2, 0), mvec.event_index(VDJ_JOI_DEL, 1, 0))
    YMIR_ASSERT2(maag.event_index(5, 1, 3, 0), mvec.event_index(VDJ_JOI_DEL, 1, 1))
    YMIR_ASSERT2(maag.event_index(5, 1, 4, 0), mvec.event_index(VDJ_JOI_DEL, 1, 2))
    YMIR_ASSERT2(maag.event_index(5, 1, 5, 0), mvec.event_index(VDJ_JOI_DEL, 1, 3))
    YMIR_ASSERT2(maag.event_index(5, 1, 6, 0), mvec.event_index(VDJ_JOI_DEL, 1, 4))
    YMIR_ASSERT2(maag.errors(2, 1, 0, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 1, 1, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 1, 2, 0), 1)
    YMIR_ASSERT2(maag.errors(2, 1, 3, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 1, 4, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 1, 5, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 1, 6, 0), 0)

    YMIR_ASSERT2(maag.event_index(5, 2, 0, 0), mvec.event_index(VDJ_JOI_DEL, 2, 0))
    YMIR_ASSERT2(maag.event_index(5, 2, 1, 0), mvec.event_index(VDJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_index(5, 2, 2, 0), mvec.event_index(VDJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(5, 2, 3, 0), mvec.event_index(VDJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(5, 2, 4, 0), mvec.event_index(VDJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(5, 2, 5, 0), mvec.event_index(VDJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(5, 2, 6, 0), mvec.event_index(VDJ_JOI_DEL, 2, 6))
    YMIR_ASSERT2(maag.errors(2, 2, 0, 0), 1)
    YMIR_ASSERT2(maag.errors(2, 2, 1, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 2, 2, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 2, 3, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 2, 4, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 2, 5, 0), 0)
    YMIR_ASSERT2(maag.errors(2, 2, 6, 0), 0)

    YMIR_ASSERT2(maag.event_index(6, 0, 0, 0), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 0, 1))
    YMIR_ASSERT2(maag.event_index(6, 0, 0, 1), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 0, 2))
    YMIR_ASSERT2(maag.event_index(6, 0, 0, 2), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 0, 0))
    YMIR_ASSERT2(maag.event_index(6, 0, 2, 0), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 2, 1))
    YMIR_ASSERT2(maag.event_index(6, 0, 2, 1), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 2, 2))
    YMIR_ASSERT2(maag.event_index(6, 0, 2, 2), mvec.event_index(VDJ_JOI_DIV_GEN, 0, 2, 0))

//    cout << "save" << endl;
//    cout << maag_builder.build(clonotype, SAVE_METADATA).fullProbability() << endl;
//    cout << "no-save" << endl;
//    cout << maag_builder.build(clonotype, NO_METADATA).fullProbability() << endl;

    // i don't want to compute by hand all this crazy matrices!
    // i've already tested chain products in previous tests!
    // ):<
    // also it's (A, B) not (A - B < eps) because this results are pretty precise on this toy example
    YMIR_ASSERT2(maag.fullProbability(0, 0, 0), maag_builder.build(clonotype, NO_METADATA, COMPUTE_ERRORS).fullProbability(0, 0, 0))  // error is here
    YMIR_ASSERT2(maag.fullProbability(1, 1, 1), maag_builder.build(clonotype, NO_METADATA, COMPUTE_ERRORS).fullProbability(1, 1, 1))
    YMIR_ASSERT2(maag.fullProbability(0, 2, 2), maag_builder.build(clonotype, NO_METADATA, COMPUTE_ERRORS).fullProbability(0, 2, 2))  // error is here

YMIR_TEST_END


YMIR_TEST_START(test_maag_builder_replace_vj)

    ModelParameterVector mvec = make_test_events_vj();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("CCCG");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCCGGG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
    seqvec2.push_back("CCGTTT");
    seqvec2.push_back("ATTT");
    seqvec2.push_back("AGGTTT");

    VDJRecombinationGenes genes("VA", alvec1, seqvec1, "JA", alvec2, seqvec2);

    MAAGBuilder maag_builder(mvec, genes);

    ClonotypeNucBuilder cl_builder;
    // CCCG.AC.GGTTT
    // poses:
    // vs: 0-1-2-3-4-5
    // js: 7-8-9-10-11
    cl_builder.setSequence("CCCGACGGTTT")
            .addVarAlignment(1, 1, 1, 4)
            .addVarAlignment(3, 1, 1, 5)
            .addJoiAlignment(1, 3, 8, 4)
            .addJoiAlignment(2, 2, 9, 3)
            .addJoiAlignment(3, 2, 7, 5);
    cl_builder.setRecombination(VJ_RECOMB);
    ClonotypeNuc clonotype = cl_builder.buildClonotype();

    MAAGnuc maag = maag_builder.build(clonotype, SAVE_METADATA, NO_ERRORS);

    ModelParameterVector mvec2 = make_test_events_vj2();

    YMIR_ASSERT(!(mvec == mvec2))

    // A .1 C .2 G .3 T .4
    // CCCGAC
    // .2 * .2 * .2 * .3 * .1 * .2 = .000048
//    cout << maag.event_probability(2, 0, 0, 0) << endl;
//    cout << mvec.event_prob(VJ_VAR_JOI_INS_LEN, 0, maag.position(6) - maag.position(0) - 1) << endl;
//    YMIR_ASSERT2(maag.event_probability(2, 0, 0, 0), mvec2.event_prob(VJ_VAR_JOI_INS_LEN, 0, maag.position(6) - maag.position(0) - 1))

    MAAGBuilder maag_builder2(mvec2, genes);

    maag_builder2.updateEventProbabilities(&maag);

    YMIR_ASSERT2(maag.nVar(), 2)
    YMIR_ASSERT2(maag.nJoi(), 3)
    YMIR_ASSERT2(maag.nDiv(), 0)

    YMIR_ASSERT2(maag.event_index(0, 0, 0, 0), mvec2.event_index(VJ_VAR_JOI_GEN, 0, 0, 0))
    YMIR_ASSERT2(maag.event_index(0, 0, 0, 1), mvec2.event_index(VJ_VAR_JOI_GEN, 0, 0, 1))
    YMIR_ASSERT2(maag.event_index(0, 0, 1, 0), mvec2.event_index(VJ_VAR_JOI_GEN, 0, 2, 0))
    YMIR_ASSERT2(maag.event_index(0, 0, 1, 2), mvec2.event_index(VJ_VAR_JOI_GEN, 0, 2, 2))

    YMIR_ASSERT2(maag.event_index(1, 0, 0, 0), mvec2.event_index(VJ_VAR_DEL, 0, 4));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 1), mvec2.event_index(VJ_VAR_DEL, 0, 3));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 2), mvec2.event_index(VJ_VAR_DEL, 0, 2));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 3), mvec2.event_index(VJ_VAR_DEL, 0, 1));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 4), mvec2.event_index(VJ_VAR_DEL, 0, 0));
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 5), 0)

    YMIR_ASSERT2(maag.event_index(1, 1, 0, 0), mvec2.event_index(VJ_VAR_DEL, 2, 6))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 1), mvec2.event_index(VJ_VAR_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 2), mvec2.event_index(VJ_VAR_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 3), mvec2.event_index(VJ_VAR_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 4), mvec2.event_index(VJ_VAR_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 5), mvec2.event_index(VJ_VAR_DEL, 2, 1))

    // A .4 C .3 G .2 T .1
    // CCCGAC
    // .3 * .3 * .3 * .2 * .4 * .3 = .000648
//    cout << maag.event_probability(2, 0, 0, 0) << endl;
//    cout << mvec2.event_prob(VJ_VAR_JOI_INS_LEN, 0, maag.position(6) - maag.position(0) - 1) << endl;
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 0), mvec2.event_index(VJ_VAR_JOI_INS_LEN, 0, maag.position(6) - maag.position(0) - 1))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 1), 0)

    YMIR_ASSERT2(maag.event_index(3, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 1, 0), mvec2.event_index(VJ_JOI_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(3, 0, 2, 0), mvec2.event_index(VJ_JOI_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(3, 0, 3, 0), mvec2.event_index(VJ_JOI_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(3, 0, 4, 0), mvec2.event_index(VJ_JOI_DEL, 0, 5))
    YMIR_ASSERT2(maag.event_index(3, 0, 5, 0), mvec2.event_index(VJ_JOI_DEL, 0, 6))

    YMIR_ASSERT2(maag.event_index(3, 2, 0, 0), mvec2.event_index(VJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_index(3, 2, 1, 0), mvec2.event_index(VJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(3, 2, 2, 0), mvec2.event_index(VJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(3, 2, 3, 0), mvec2.event_index(VJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(3, 2, 4, 0), mvec2.event_index(VJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(3, 2, 5, 0), mvec2.event_index(VJ_JOI_DEL, 2, 6))

YMIR_TEST_END


YMIR_TEST_START(test_maag_builder_replace_vdj)

    ModelParameterVector mvec = make_test_events_vdj();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("CCCG");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCCGGG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
    seqvec2.push_back("CCGTTT");
    seqvec2.push_back("ATTT");
    seqvec2.push_back("AGGTTT");

    vector<string> alvec3;
    vector<string> seqvec3;
    alvec3.push_back("Dseg1");
    alvec3.push_back("Dseg2");
    alvec3.push_back("Dseg3");
    seqvec3.push_back("GTTT");
    seqvec3.push_back("ACCGGT");
    seqvec3.push_back("CCCGGAC");

    VDJRecombinationGenes genes("VB", alvec1, seqvec1, "JB", alvec2, seqvec2, "DB", alvec3, seqvec3);

    MAAGBuilder maag_builder(mvec, genes);

    ClonotypeNucBuilder cl_builder;
    /*
     D1:
       CCCGACGGTTT
             .GTTT
             .GTT.T
             G.TTT
     D2:
       CCCGACGGTTT
      A.CCG.GT
         AC.CGGT

     D3:
       CCCGACGGTTT
     CCCG.GAC
       CCCG.GAC
         CC.CGG.AC
    */
    cl_builder.setSequence("CCCGACGGTTT")
            .addVarAlignment(1, 1, 1, 4)
            .addVarAlignment(3, 1, 1, 5)
            .addJoiAlignment(1, 3, 8, 4)
            .addJoiAlignment(2, 2, 9, 3)
            .addJoiAlignment(3, 2, 7, 5)
            .addDivAlignment(2, 2, 2, 3)
            .addDivAlignment(2, 3, 6, 4)
            .addDivAlignment(3, 5, 4, 3)
            .addDivAlignment(3, 1, 1, 4)
            .addDivAlignment(3, 3, 6, 3)
            .addDivAlignment(1, 1, 8, 4);
    cl_builder.setRecombination(VDJ_RECOMB);
    ClonotypeNuc clonotype = cl_builder.buildClonotype();

    MAAGnuc maag = maag_builder.build(clonotype, SAVE_METADATA, NO_ERRORS);

    ModelParameterVector mvec2 = make_test_events_vdj2();

    YMIR_ASSERT(!(mvec == mvec2))

    MAAGBuilder maag_builder2(mvec2, genes);
    maag_builder2.updateEventProbabilities(&maag);

    YMIR_ASSERT2(maag.nVar(), 2)
    YMIR_ASSERT2(maag.nJoi(), 3)
    YMIR_ASSERT2(maag.nDiv(), 3)

    YMIR_ASSERT2(maag.event_index(0, 0, 0, 0), mvec2.event_index(VDJ_VAR_GEN, 0, 0))
    YMIR_ASSERT2(maag.event_index(0, 1, 0, 0), mvec2.event_index(VDJ_VAR_GEN, 0, 2))

    YMIR_ASSERT2(maag.event_index(1, 0, 0, 0), mvec2.event_index(VDJ_VAR_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 1), mvec2.event_index(VDJ_VAR_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 2), mvec2.event_index(VDJ_VAR_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 3), mvec2.event_index(VDJ_VAR_DEL, 0, 1))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 4), mvec2.event_index(VDJ_VAR_DEL, 0, 0))
    YMIR_ASSERT2(maag.event_index(1, 0, 0, 5), 0)

    YMIR_ASSERT2(maag.event_index(1, 1, 0, 0), mvec2.event_index(VDJ_VAR_DEL, 2, 6))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 1), mvec2.event_index(VDJ_VAR_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 2), mvec2.event_index(VDJ_VAR_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 3), mvec2.event_index(VDJ_VAR_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 4), mvec2.event_index(VDJ_VAR_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(1, 1, 0, 5), mvec2.event_index(VDJ_VAR_DEL, 2, 1))

    // seq: 1 2 4 6 7 8 9
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 0), mvec2.event_index(VDJ_VAR_DIV_INS_LEN, 0, 0))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 1), mvec2.event_index(VDJ_VAR_DIV_INS_LEN, 0, 1))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 2), mvec2.event_index(VDJ_VAR_DIV_INS_LEN, 0, 3))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 3), mvec2.event_index(VDJ_VAR_DIV_INS_LEN, 0, 5))
    YMIR_ASSERT2(maag.event_index(2, 0, 2, 2), mvec2.event_index(VDJ_VAR_DIV_INS_LEN, 0, 1))

    // rows 1 2 3 4 5 | cols 3 4 5 6 7
    // TODO: check D deletions indices
    YMIR_ASSERT2(maag.event_index(3, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 0, 1), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 1, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 1, 1), mvec.event_index(VDJ_DIV_DEL, 1, 2, 2))
    YMIR_ASSERT2(maag.event_index(3, 0, 2, 2), 0)
    YMIR_ASSERT2(maag.event_index(3, 0, 3, 3), mvec.event_index(VDJ_DIV_DEL, 1, 3, 1))
    YMIR_ASSERT2(maag.event_index(3, 0, 3, 4), mvec.event_index(VDJ_DIV_DEL, 1, 3, 0))
    YMIR_ASSERT2(maag.event_index(3, 0, 4, 4), mvec.event_index(VDJ_DIV_DEL, 1, 4, 0))

    YMIR_ASSERT2(maag.event_index(3, 2, 5, 5), mvec.event_index(VDJ_DIV_DEL, 0, 0, 3))
    YMIR_ASSERT2(maag.event_index(3, 2, 6, 5), 0)
    YMIR_ASSERT2(maag.event_index(3, 2, 5, 6), mvec.event_index(VDJ_DIV_DEL, 0, 1, 0))

    // row: 3 4 6 7 8 10 11
    // col: 7 8 9 10 11 [12]
    YMIR_ASSERT2(maag.event_index(4, 0, 0, 0), mvec2.event_index(VDJ_DIV_JOI_INS_LEN, 0, 3))
    YMIR_ASSERT2(maag.event_index(4, 0, 1, 1), mvec2.event_index(VDJ_DIV_JOI_INS_LEN, 0, 3))
    YMIR_ASSERT2(maag.event_index(4, 0, 3, 0), 0)

    YMIR_ASSERT2(maag.event_index(5, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(5, 0, 1, 0), mvec2.event_index(VDJ_JOI_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(5, 0, 2, 0), mvec2.event_index(VDJ_JOI_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(5, 0, 3, 0), mvec2.event_index(VDJ_JOI_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(5, 0, 4, 0), mvec2.event_index(VDJ_JOI_DEL, 0, 5))
    YMIR_ASSERT2(maag.event_index(5, 0, 5, 0), mvec2.event_index(VDJ_JOI_DEL, 0, 6))

    YMIR_ASSERT2(maag.event_index(5, 2, 0, 0), mvec2.event_index(VDJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_index(5, 2, 1, 0), mvec2.event_index(VDJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(5, 2, 2, 0), mvec2.event_index(VDJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(5, 2, 3, 0), mvec2.event_index(VDJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(5, 2, 4, 0), mvec2.event_index(VDJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(5, 2, 5, 0), mvec2.event_index(VDJ_JOI_DEL, 2, 6))

    YMIR_ASSERT2(maag.event_index(6, 0, 0, 0), mvec2.event_index(VDJ_JOI_DIV_GEN, 0, 0, 1))
    YMIR_ASSERT2(maag.event_index(6, 0, 0, 1), mvec2.event_index(VDJ_JOI_DIV_GEN, 0, 0, 2))
    YMIR_ASSERT2(maag.event_index(6, 0, 0, 2), mvec2.event_index(VDJ_JOI_DIV_GEN, 0, 0, 0))
    YMIR_ASSERT2(maag.event_index(6, 0, 2, 0), mvec2.event_index(VDJ_JOI_DIV_GEN, 0, 2, 1))
    YMIR_ASSERT2(maag.event_index(6, 0, 2, 1), mvec2.event_index(VDJ_JOI_DIV_GEN, 0, 2, 2))
    YMIR_ASSERT2(maag.event_index(6, 0, 2, 2), mvec2.event_index(VDJ_JOI_DIV_GEN, 0, 2, 0))

YMIR_TEST_END


YMIR_TEST_START(test_maag_vj_aa)

    ModelParameterVector mvec = make_test_events_vj3();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1"); seqvec1.push_back("TGTGC");
    alvec1.push_back("Vseg2"); seqvec1.push_back("TGCG");
    alvec1.push_back("Vseg3"); seqvec1.push_back("TGA");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1"); seqvec2.push_back("TTT");
    alvec2.push_back("Jseg2"); seqvec2.push_back("ATTC");
    alvec2.push_back("Jseg3"); seqvec2.push_back("AATT");

    VDJRecombinationGenes genes("VA", alvec1, seqvec1, "JA", alvec2, seqvec2);

    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)] = .1;
    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 1)] = .2;
    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 2)] = .3;
    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 3)] = .4;

    MAAGBuilder maag_builder(mvec, genes);

    ClonotypeNucBuilder cl_builder;
    // C: TGT TGC
    // A: GCT GCC GCA GCG
    // S: TCT TCC TCA TCG AGT AGC
    // F: TTT TTC
    CDR3AminoAcidAligner aligner(genes, VDJAlignerParameters(3,
                                                             VDJAlignmentEventScore(AlignmentEventScore(1, -1, 1),
                                                                                    AlignmentEventScore(1, -1, 1),
                                                                                    AlignmentEventScore(1, -1, 1)),
                                                             VDJAlignmentScoreThreshold(1, 2, 1))); // VDJAlignmentScoreThreshold(2, 2, 2))); - error with this string
    aligner.setSequence("CASF").setRecombination(VJ_RECOMB);
    YMIR_ASSERT(aligner.alignVar())
    YMIR_ASSERT(aligner.alignJoi())

    ClonotypeAA clonotype = aligner.buildClonotype();

    YMIR_ASSERT2(clonotype.getVarLen(0), 5)
    YMIR_ASSERT2(clonotype.getVarLen(1), 4)
    YMIR_ASSERT2(clonotype.getVarLen(2), 2)
    YMIR_ASSERT2(clonotype.getVarCodon(0, 1), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERT2(clonotype.getVarCodon(0, 2), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERT2(clonotype.getVarCodon(0, 3), compute_codon_hash({true, false, false, false, false, false}, 0))
    YMIR_ASSERT2(clonotype.getVarCodon(0, 4), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERT2(clonotype.getVarCodon(0, 5), compute_codon_hash({true, true, true, true, false, false}, 0))

    YMIR_ASSERT2(clonotype.getJoiLen(0), 3)
    YMIR_ASSERT2(clonotype.getJoiLen(1), 4)
    YMIR_ASSERT2(clonotype.getJoiLen(2), 2)
    YMIR_ASSERT2(clonotype.getJoiSeqStart(0), 10)
    YMIR_ASSERT2(clonotype.getJoiSeqEnd(0), 12)
    YMIR_ASSERT2(clonotype.getJoiGeneStart(0), 1)
    YMIR_ASSERT2(clonotype.getJoiGeneEnd(0), 3)
    YMIR_ASSERT2(clonotype.getJoiCodon(0, 1), compute_codon_hash({true, false, false, false, false, false}, 0))
    YMIR_ASSERT2(clonotype.getJoiCodon(0, 2), compute_codon_hash({true, false, false, false, false, false}, 0))
    YMIR_ASSERT2(clonotype.getJoiCodon(0, 3), compute_codon_hash({true, false, false, false, false, false}, 0))

    MAAGaa maag = maag_builder.build(clonotype);

    YMIR_ASSERT2(maag.nVar(), 3)
    YMIR_ASSERT2(maag.nDiv(), 0)
    YMIR_ASSERT2(maag.nJoi(), 3)

    YMIR_ASSERT2(maag.n_poses(), 11)

    YMIR_ASSERT2(maag.position(0), 0)
    YMIR_ASSERT2(maag.position(1), 1)
    YMIR_ASSERT2(maag.position(2), 2)
    YMIR_ASSERT2(maag.position(3), 3)
    YMIR_ASSERT2(maag.position(4), 4)
    YMIR_ASSERT2(maag.position(5), 5)

    YMIR_ASSERT2(maag.position(6), 9)
    YMIR_ASSERT2(maag.position(7), 10)
    YMIR_ASSERT2(maag.position(8), 11)
    YMIR_ASSERT2(maag.position(9), 12)
    YMIR_ASSERT2(maag.position(10), 13)

    std::vector<string> rev_nuc =  {"TGCGCAAGCTTC", "TGTGCAAGCTTC", "TGCGCCAGCTTC", "TGTGCCAGCTTC", "TGCGCGAGCTTC", "TGTGCGAGCTTC",
                                    "TGCGCTAGCTTC", "TGTGCTAGCTTC", "TGCGCAAGTTTC", "TGTGCAAGTTTC", "TGCGCCAGTTTC", "TGTGCCAGTTTC",
                                    "TGCGCGAGTTTC", "TGTGCGAGTTTC", "TGCGCTAGTTTC", "TGTGCTAGTTTC", "TGCGCATCATTC", "TGTGCATCATTC",
                                    "TGCGCCTCATTC", "TGTGCCTCATTC", "TGCGCGTCATTC", "TGTGCGTCATTC", "TGCGCTTCATTC", "TGTGCTTCATTC",
                                    "TGCGCATCCTTC", "TGTGCATCCTTC", "TGCGCCTCCTTC", "TGTGCCTCCTTC", "TGCGCGTCCTTC", "TGTGCGTCCTTC",
                                    "TGCGCTTCCTTC", "TGTGCTTCCTTC", "TGCGCATCGTTC", "TGTGCATCGTTC", "TGCGCCTCGTTC", "TGTGCCTCGTTC",
                                    "TGCGCGTCGTTC", "TGTGCGTCGTTC", "TGCGCTTCGTTC", "TGTGCTTCGTTC", "TGCGCATCTTTC", "TGTGCATCTTTC",
                                    "TGCGCCTCTTTC", "TGTGCCTCTTTC", "TGCGCGTCTTTC", "TGTGCGTCTTTC", "TGCGCTTCTTTC", "TGTGCTTCTTTC",
                                    "TGCGCAAGCTTT", "TGTGCAAGCTTT", "TGCGCCAGCTTT", "TGTGCCAGCTTT", "TGCGCGAGCTTT", "TGTGCGAGCTTT",
                                    "TGCGCTAGCTTT", "TGTGCTAGCTTT", "TGCGCAAGTTTT", "TGTGCAAGTTTT", "TGCGCCAGTTTT", "TGTGCCAGTTTT",
                                    "TGCGCGAGTTTT", "TGTGCGAGTTTT", "TGCGCTAGTTTT", "TGTGCTAGTTTT", "TGCGCATCATTT", "TGTGCATCATTT",
                                    "TGCGCCTCATTT", "TGTGCCTCATTT", "TGCGCGTCATTT", "TGTGCGTCATTT", "TGCGCTTCATTT", "TGTGCTTCATTT",
                                    "TGCGCATCCTTT", "TGTGCATCCTTT", "TGCGCCTCCTTT", "TGTGCCTCCTTT", "TGCGCGTCCTTT", "TGTGCGTCCTTT",
                                    "TGCGCTTCCTTT", "TGTGCTTCCTTT", "TGCGCATCGTTT", "TGTGCATCGTTT", "TGCGCCTCGTTT", "TGTGCCTCGTTT",
                                    "TGCGCGTCGTTT", "TGTGCGTCGTTT", "TGCGCTTCGTTT", "TGTGCTTCGTTT", "TGCGCATCTTTT", "TGTGCATCTTTT",
                                    "TGCGCCTCTTTT", "TGTGCCTCTTTT", "TGCGCGTCTTTT", "TGTGCGTCTTTT", "TGCGCTTCTTTT", "TGTGCTTCTTTT"};

    std::vector<ClonotypeNuc> clonotype_vec;
    NaiveCDR3NucleotideAligner naligner(genes, VDJAlignerParameters(3));
    for (auto &seq: rev_nuc) {
        naligner.setSequence(seq).setRecombination(VJ_RECOMB);
        YMIR_ASSERT(naligner.alignVar())
        YMIR_ASSERT(naligner.alignJoi())
        clonotype_vec.push_back(naligner.buildClonotype());
    }

    prob_t sum_prob = 0;
    std::array<prob_t, 9> prob_vec = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < clonotype_vec.size(); ++i) {
        sum_prob += maag_builder.buildAndCompute(clonotype_vec[i], NO_ERRORS);
        for (int k = 0; k < 3; ++k) {
            for (int j = 0; j < 3; ++j) {
                prob_vec[3*k + j] += maag_builder.build(clonotype_vec[i], NO_METADATA, NO_ERRORS).fullProbability(k, j);
            }
        }
    }

    YMIR_ASSERT3(maag.fullProbability(0, 0), prob_vec[0])
    YMIR_ASSERT3(maag.fullProbability(0, 1), prob_vec[1])
    YMIR_ASSERT3(maag.fullProbability(0, 2), prob_vec[2])
    YMIR_ASSERT3(maag.fullProbability(1, 0), prob_vec[3])
    YMIR_ASSERT3(maag.fullProbability(1, 1), prob_vec[4])
    YMIR_ASSERT3(maag.fullProbability(1, 2), prob_vec[5])
    YMIR_ASSERT3(maag.fullProbability(2, 0), prob_vec[6])
    YMIR_ASSERT3(maag.fullProbability(2, 1), prob_vec[7])
    YMIR_ASSERT3(maag.fullProbability(2, 2), prob_vec[8])

    YMIR_ASSERT3(sum_prob, prob_vec[0] + prob_vec[1] + prob_vec[2] + prob_vec[3] + prob_vec[4] + prob_vec[5] + prob_vec[6] + prob_vec[7] + prob_vec[8])
    YMIR_ASSERT3(maag.fullProbability(), prob_vec[0] + prob_vec[1] + prob_vec[2] + prob_vec[3] + prob_vec[4] + prob_vec[5] + prob_vec[6] + prob_vec[7] + prob_vec[8])
    YMIR_ASSERT3(maag.fullProbability(), sum_prob)
    YMIR_ASSERT3(maag.fullProbability(), maag_builder.buildAndCompute(clonotype))

YMIR_TEST_END


YMIR_TEST_START(test_maag_vdj_aa_very_simple)

//    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vdj_model3/");

    ModelParameterVector mvec = make_test_events_vdj4();

    for (int i = 0; i < mvec.size(); ++i) {
        mvec[i] = 1;
    }
//    mvec.normaliseEventFamilies();

    vector<string> alvec1;
    vector<string> seqvec1;
//    alvec1.push_back("Vseg1"); seqvec1.push_back("AAAAA");
    alvec1.push_back("Vseg1"); seqvec1.push_back("AATAA");
    alvec1.push_back("Vseg2"); seqvec1.push_back("GGGG");
    alvec1.push_back("Vseg3"); seqvec1.push_back("GGG");

    vector<string> alvec2;
    vector<string> seqvec2;
//    alvec2.push_back("Jseg1"); seqvec2.push_back("CTT");
//    alvec2.push_back("Jseg1"); seqvec2.push_back("TTT");
    alvec2.push_back("Jseg1"); seqvec2.push_back("TGG");
    alvec2.push_back("Jseg2"); seqvec2.push_back("GGGG");
    alvec2.push_back("Jseg3"); seqvec2.push_back("GGGG");

    vector<string> alvec3;
    vector<string> seqvec3;
//    alvec3.push_back("Dseg1"); seqvec3.push_back("CCC");
//    alvec3.push_back("Dseg1"); seqvec3.push_back("ACT");
    alvec3.push_back("Dseg1"); seqvec3.push_back("GGG");
    alvec3.push_back("Dseg2"); seqvec3.push_back("AAA");

    VDJRecombinationGenes genes("VB", alvec1, seqvec1, "JB", alvec2, seqvec2, "DB", alvec3, seqvec3);

    MAAGBuilder maag_builder(mvec, genes);


    auto aligner_params = VDJAlignerParameters(3,
                                               VDJAlignmentEventScore(AlignmentEventScore(1, -1, 1),
                                                                      AlignmentEventScore(1, -1, 1),
                                                                      AlignmentEventScore(1, -1, 1)),
                                               VDJAlignmentScoreThreshold(3, 2, 3));
    ClonotypeNucBuilder cl_builder;
    CDR3AminoAcidAligner aligner(genes, aligner_params);
//    aligner.setSequence("KPF").setRecombination(VDJ_RECOMB);
//    aligner.setSequence("KL").setRecombination(VDJ_RECOMB);
    aligner.setSequence("NGW").setRecombination(VDJ_RECOMB);
    YMIR_ASSERT(aligner.alignVar())
    YMIR_ASSERT(aligner.alignDiv())
    YMIR_ASSERT(aligner.alignJoi())

//    std::vector<string> rev_nuc = {"AAACCATTC", "AAGCCATTC", "AAACCCTTC", "AAGCCCTTC", "AAACCGTTC", "AAGCCGTTC", "AAACCTTTC", "AAGCCTTTC", "AAACCATTT", "AAGCCATTT", "AAACCCTTT", "AAGCCCTTT", "AAACCGTTT", "AAGCCGTTT", "AAACCTTTT", "AAGCCTTTT"};

//    std::vector<string> rev_nuc = {"AAACTA", "AAGCTA", "AAACTC", "AAGCTC", "AAACTG", "AAGCTG", "AAACTT", "AAGCTT", "AAATTA", "AAGTTA", "AAATTG", "AAGTTG"};

    std::vector<string> rev_nuc = {"AATGGATGG", "AATGGCTGG", "AATGGGTGG", "AATGGTTGG", "AACGGATGG", "AACGGCTGG", "AACGGGTGG", "AACGGTTGG"};

//    std::vector<string> rev_nuc = {"AAACCCTTT", "AAGCCCTTT"};

    std::cout << "aa building" << std::endl;
    ClonotypeAA clonotype = aligner.buildClonotype();

    YMIR_ASSERT2((int) clonotype.nVar(), 1)
    YMIR_ASSERT2((int) clonotype.nDiv(), 1)
    YMIR_ASSERT2((int) clonotype.nJoi(), 1)

    MAAGaa maag = maag_builder.build(clonotype);
    maag.fullProbability();
    maag.printCodons();
    maag.print();

    std::vector<ClonotypeNuc> clonotype_vec;
    NaiveCDR3NucleotideAligner naligner(genes, aligner_params);

    std::cout << "nuc building" << std::endl;
    for (auto &seq: rev_nuc) {
        naligner.clear();
        naligner.setSequence(seq).setRecombination(VDJ_RECOMB);
        if (naligner.alignVar() && naligner.alignDiv() && naligner.alignJoi()) {
            clonotype_vec.push_back(naligner.buildClonotype());
        }
    }

    prob_t sum_prob = 0;
    std::vector<prob_t> prob_vec;
    for (int i = 0; i < clonotype_vec.size(); ++i) {
        sum_prob += maag_builder.buildAndCompute(clonotype_vec[i], NO_ERRORS);
        auto tmp_maag = maag_builder.build(clonotype_vec[i], NO_METADATA, NO_ERRORS);
        for (int v_i = 0; v_i < tmp_maag.nVar(); ++v_i) {
            for (int d_i = 0; d_i < tmp_maag.nDiv(); ++d_i) {
                for (int j_i = 0; j_i < tmp_maag.nJoi(); ++j_i) {
                    prob_vec.push_back(tmp_maag.fullProbability(v_i, d_i, j_i));
                }
                tmp_maag.print();
            }
        }
    }

    std::cout << "size:" << (int) prob_vec.size() << std::endl;
    YMIR_ASSERT3(sum_prob, std::accumulate(prob_vec.begin(), prob_vec.end(), .0))
    YMIR_ASSERT3(maag.fullProbability(), std::accumulate(prob_vec.begin(), prob_vec.end(), .0))
    YMIR_ASSERT3(maag.fullProbability(), sum_prob)
    YMIR_ASSERT3(maag.fullProbability(), maag_builder.buildAndCompute(clonotype))

YMIR_TEST_END


YMIR_TEST_START(test_maag_vdj_aa_simple_D)

    ModelParameterVector mvec = make_test_events_vdj3();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1"); seqvec1.push_back("TGTGC");
    alvec1.push_back("Vseg2"); seqvec1.push_back("TGCG");
    alvec1.push_back("Vseg3"); seqvec1.push_back("TGA");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1"); seqvec2.push_back("TTT");
    alvec2.push_back("Jseg2"); seqvec2.push_back("ATTC");
    alvec2.push_back("Jseg3"); seqvec2.push_back("AATT");

    vector<string> alvec3;
    vector<string> seqvec3;
    alvec3.push_back("Dseg1"); seqvec3.push_back("AATT");
    alvec3.push_back("Dseg2"); seqvec3.push_back("GGGG");

    VDJRecombinationGenes genes("VB", alvec1, seqvec1, "JB", alvec2, seqvec2, "DB", alvec3, seqvec3);

    MAAGBuilder maag_builder(mvec, genes);


    ClonotypeNucBuilder cl_builder;
    CDR3AminoAcidAligner aligner(genes, VDJAlignerParameters(3,
                                                             VDJAlignmentEventScore(AlignmentEventScore(1, -1, 1),
                                                                                    AlignmentEventScore(1, -1, 1),
                                                                                    AlignmentEventScore(1, -1, 1)),
                                                             VDJAlignmentScoreThreshold(2, 2, 2)));
    aligner.setSequence("CAQF").setRecombination(VDJ_RECOMB);
    YMIR_ASSERT(aligner.alignVar())
    YMIR_ASSERT(aligner.alignDiv())
    YMIR_ASSERT(aligner.alignJoi())

    std::vector<string> rev_nuc = {"TGCGCACAATTC", "TGTGCACAATTC", "TGCGCCCAATTC", "TGTGCCCAATTC", "TGCGCGCAATTC", "TGTGCGCAATTC", "TGCGCTCAATTC", "TGTGCTCAATTC", "TGCGCACAGTTC", "TGTGCACAGTTC", "TGCGCCCAGTTC", "TGTGCCCAGTTC", "TGCGCGCAGTTC", "TGTGCGCAGTTC", "TGCGCTCAGTTC", "TGTGCTCAGTTC", "TGCGCACAATTT", "TGTGCACAATTT", "TGCGCCCAATTT", "TGTGCCCAATTT", "TGCGCGCAATTT", "TGTGCGCAATTT", "TGCGCTCAATTT", "TGTGCTCAATTT", "TGCGCACAGTTT", "TGTGCACAGTTT", "TGCGCCCAGTTT", "TGTGCCCAGTTT", "TGCGCGCAGTTT", "TGTGCGCAGTTT", "TGCGCTCAGTTT", "TGTGCTCAGTTT"};

    ClonotypeAA clonotype = aligner.buildClonotype();

    YMIR_ASSERT2(clonotype.nVar(), 3)
    YMIR_ASSERT2(clonotype.nDiv(), 1)
    YMIR_ASSERT2(clonotype.nJoi(), 3)

    MAAGaa maag = maag_builder.build(clonotype);

    std::vector<ClonotypeNuc> clonotype_vec;
    NaiveCDR3NucleotideAligner naligner(genes, VDJAlignerParameters(3,
                                                                    VDJAlignmentEventScore(AlignmentEventScore(1, -1, 1),
                                                                                           AlignmentEventScore(1, -1, 1),
                                                                                           AlignmentEventScore(1, -1, 1)),
                                                                    VDJAlignmentScoreThreshold(2, 2, 2)));
    for (auto &seq: rev_nuc) {
        naligner.clear();
        naligner.setSequence(seq).setRecombination(VDJ_RECOMB);
        if (naligner.alignVar() && naligner.alignDiv() && naligner.alignJoi()) {
            clonotype_vec.push_back(naligner.buildClonotype());
        }
    }

    prob_t sum_prob = 0;
    std::vector<prob_t> prob_vec;
    for (int i = 0; i < clonotype_vec.size(); ++i) {
        sum_prob += maag_builder.buildAndCompute(clonotype_vec[i], NO_ERRORS);
        auto tmp_maag = maag_builder.build(clonotype_vec[i], NO_METADATA, NO_ERRORS);
        for (int v_i = 0; v_i < tmp_maag.nVar(); ++v_i) {
            for (int d_i = 0; d_i < tmp_maag.nDiv(); ++d_i) {
                for (int j_i = 0; j_i < tmp_maag.nJoi(); ++j_i) {
                    prob_vec.push_back(tmp_maag.fullProbability(v_i, d_i, j_i));
                }
            }
        }
    }

    // result: 2.11458e-06
    YMIR_ASSERT3(sum_prob, std::accumulate(prob_vec.begin(), prob_vec.end(), .0))
    YMIR_ASSERT3(maag.fullProbability(), std::accumulate(prob_vec.begin(), prob_vec.end(), .0))
    YMIR_ASSERT3(maag.fullProbability(), sum_prob)
    YMIR_ASSERT3(maag.fullProbability(), maag_builder.buildAndCompute(clonotype))

YMIR_TEST_END


YMIR_TEST_START(test_maag_vdj_aa)

    ModelParameterVector mvec = make_test_events_vdj3();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1"); seqvec1.push_back("TGTGC");
    alvec1.push_back("Vseg2"); seqvec1.push_back("TGCG");
    alvec1.push_back("Vseg3"); seqvec1.push_back("TGA");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1"); seqvec2.push_back("TTT");
    alvec2.push_back("Jseg2"); seqvec2.push_back("ATTC");
    alvec2.push_back("Jseg3"); seqvec2.push_back("AATT");

    vector<string> alvec3;
    vector<string> seqvec3;
    alvec3.push_back("Dseg1"); seqvec3.push_back("GGCA");
    alvec3.push_back("Dseg2"); seqvec3.push_back("GGTT");

    VDJRecombinationGenes genes("VB", alvec1, seqvec1, "JB", alvec2, seqvec2, "DB", alvec3, seqvec3);

    MAAGBuilder maag_builder(mvec, genes);


    ClonotypeNucBuilder cl_builder;
    CDR3AminoAcidAligner aligner(genes, VDJAlignerParameters(3,
                                                             VDJAlignmentEventScore(AlignmentEventScore(1, -1, 1),
                                                                                    AlignmentEventScore(1, -1, 1),
                                                                                    AlignmentEventScore(1, -1, 1)),
                                                             VDJAlignmentScoreThreshold(2, 2, 2)));
    aligner.setSequence("CAGVQF").setRecombination(VDJ_RECOMB);
    YMIR_ASSERT(aligner.alignVar())
    YMIR_ASSERT(aligner.alignDiv())
    YMIR_ASSERT(aligner.alignJoi())

    std::vector<string> rev_nuc = {"TGCGCAGGAGTACAATTC", "TGTGCAGGAGTACAATTC", "TGCGCCGGAGTACAATTC", "TGTGCCGGAGTACAATTC", "TGCGCGGGAGTACAATTC", "TGTGCGGGAGTACAATTC", "TGCGCTGGAGTACAATTC", "TGTGCTGGAGTACAATTC", "TGCGCAGGCGTACAATTC", "TGTGCAGGCGTACAATTC", "TGCGCCGGCGTACAATTC", "TGTGCCGGCGTACAATTC", "TGCGCGGGCGTACAATTC", "TGTGCGGGCGTACAATTC", "TGCGCTGGCGTACAATTC", "TGTGCTGGCGTACAATTC", "TGCGCAGGGGTACAATTC", "TGTGCAGGGGTACAATTC", "TGCGCCGGGGTACAATTC", "TGTGCCGGGGTACAATTC", "TGCGCGGGGGTACAATTC", "TGTGCGGGGGTACAATTC", "TGCGCTGGGGTACAATTC", "TGTGCTGGGGTACAATTC", "TGCGCAGGTGTACAATTC", "TGTGCAGGTGTACAATTC", "TGCGCCGGTGTACAATTC", "TGTGCCGGTGTACAATTC", "TGCGCGGGTGTACAATTC", "TGTGCGGGTGTACAATTC", "TGCGCTGGTGTACAATTC", "TGTGCTGGTGTACAATTC", "TGCGCAGGAGTCCAATTC", "TGTGCAGGAGTCCAATTC", "TGCGCCGGAGTCCAATTC", "TGTGCCGGAGTCCAATTC", "TGCGCGGGAGTCCAATTC", "TGTGCGGGAGTCCAATTC", "TGCGCTGGAGTCCAATTC", "TGTGCTGGAGTCCAATTC", "TGCGCAGGCGTCCAATTC", "TGTGCAGGCGTCCAATTC", "TGCGCCGGCGTCCAATTC", "TGTGCCGGCGTCCAATTC", "TGCGCGGGCGTCCAATTC", "TGTGCGGGCGTCCAATTC", "TGCGCTGGCGTCCAATTC", "TGTGCTGGCGTCCAATTC", "TGCGCAGGGGTCCAATTC", "TGTGCAGGGGTCCAATTC", "TGCGCCGGGGTCCAATTC", "TGTGCCGGGGTCCAATTC", "TGCGCGGGGGTCCAATTC", "TGTGCGGGGGTCCAATTC", "TGCGCTGGGGTCCAATTC", "TGTGCTGGGGTCCAATTC", "TGCGCAGGTGTCCAATTC", "TGTGCAGGTGTCCAATTC", "TGCGCCGGTGTCCAATTC", "TGTGCCGGTGTCCAATTC", "TGCGCGGGTGTCCAATTC", "TGTGCGGGTGTCCAATTC", "TGCGCTGGTGTCCAATTC", "TGTGCTGGTGTCCAATTC", "TGCGCAGGAGTGCAATTC", "TGTGCAGGAGTGCAATTC", "TGCGCCGGAGTGCAATTC", "TGTGCCGGAGTGCAATTC", "TGCGCGGGAGTGCAATTC", "TGTGCGGGAGTGCAATTC", "TGCGCTGGAGTGCAATTC", "TGTGCTGGAGTGCAATTC", "TGCGCAGGCGTGCAATTC", "TGTGCAGGCGTGCAATTC", "TGCGCCGGCGTGCAATTC", "TGTGCCGGCGTGCAATTC", "TGCGCGGGCGTGCAATTC", "TGTGCGGGCGTGCAATTC", "TGCGCTGGCGTGCAATTC", "TGTGCTGGCGTGCAATTC", "TGCGCAGGGGTGCAATTC", "TGTGCAGGGGTGCAATTC", "TGCGCCGGGGTGCAATTC", "TGTGCCGGGGTGCAATTC", "TGCGCGGGGGTGCAATTC", "TGTGCGGGGGTGCAATTC", "TGCGCTGGGGTGCAATTC", "TGTGCTGGGGTGCAATTC", "TGCGCAGGTGTGCAATTC", "TGTGCAGGTGTGCAATTC", "TGCGCCGGTGTGCAATTC", "TGTGCCGGTGTGCAATTC", "TGCGCGGGTGTGCAATTC", "TGTGCGGGTGTGCAATTC", "TGCGCTGGTGTGCAATTC", "TGTGCTGGTGTGCAATTC", "TGCGCAGGAGTTCAATTC", "TGTGCAGGAGTTCAATTC", "TGCGCCGGAGTTCAATTC", "TGTGCCGGAGTTCAATTC", "TGCGCGGGAGTTCAATTC", "TGTGCGGGAGTTCAATTC", "TGCGCTGGAGTTCAATTC", "TGTGCTGGAGTTCAATTC", "TGCGCAGGCGTTCAATTC", "TGTGCAGGCGTTCAATTC", "TGCGCCGGCGTTCAATTC", "TGTGCCGGCGTTCAATTC", "TGCGCGGGCGTTCAATTC", "TGTGCGGGCGTTCAATTC", "TGCGCTGGCGTTCAATTC", "TGTGCTGGCGTTCAATTC", "TGCGCAGGGGTTCAATTC", "TGTGCAGGGGTTCAATTC", "TGCGCCGGGGTTCAATTC", "TGTGCCGGGGTTCAATTC", "TGCGCGGGGGTTCAATTC", "TGTGCGGGGGTTCAATTC", "TGCGCTGGGGTTCAATTC", "TGTGCTGGGGTTCAATTC", "TGCGCAGGTGTTCAATTC", "TGTGCAGGTGTTCAATTC", "TGCGCCGGTGTTCAATTC", "TGTGCCGGTGTTCAATTC", "TGCGCGGGTGTTCAATTC", "TGTGCGGGTGTTCAATTC", "TGCGCTGGTGTTCAATTC", "TGTGCTGGTGTTCAATTC", "TGCGCAGGAGTACAGTTC", "TGTGCAGGAGTACAGTTC", "TGCGCCGGAGTACAGTTC", "TGTGCCGGAGTACAGTTC", "TGCGCGGGAGTACAGTTC", "TGTGCGGGAGTACAGTTC", "TGCGCTGGAGTACAGTTC", "TGTGCTGGAGTACAGTTC", "TGCGCAGGCGTACAGTTC", "TGTGCAGGCGTACAGTTC", "TGCGCCGGCGTACAGTTC", "TGTGCCGGCGTACAGTTC", "TGCGCGGGCGTACAGTTC", "TGTGCGGGCGTACAGTTC", "TGCGCTGGCGTACAGTTC", "TGTGCTGGCGTACAGTTC", "TGCGCAGGGGTACAGTTC", "TGTGCAGGGGTACAGTTC", "TGCGCCGGGGTACAGTTC", "TGTGCCGGGGTACAGTTC", "TGCGCGGGGGTACAGTTC", "TGTGCGGGGGTACAGTTC", "TGCGCTGGGGTACAGTTC", "TGTGCTGGGGTACAGTTC", "TGCGCAGGTGTACAGTTC", "TGTGCAGGTGTACAGTTC", "TGCGCCGGTGTACAGTTC", "TGTGCCGGTGTACAGTTC", "TGCGCGGGTGTACAGTTC", "TGTGCGGGTGTACAGTTC", "TGCGCTGGTGTACAGTTC", "TGTGCTGGTGTACAGTTC", "TGCGCAGGAGTCCAGTTC", "TGTGCAGGAGTCCAGTTC", "TGCGCCGGAGTCCAGTTC", "TGTGCCGGAGTCCAGTTC", "TGCGCGGGAGTCCAGTTC", "TGTGCGGGAGTCCAGTTC", "TGCGCTGGAGTCCAGTTC", "TGTGCTGGAGTCCAGTTC", "TGCGCAGGCGTCCAGTTC", "TGTGCAGGCGTCCAGTTC", "TGCGCCGGCGTCCAGTTC", "TGTGCCGGCGTCCAGTTC", "TGCGCGGGCGTCCAGTTC", "TGTGCGGGCGTCCAGTTC", "TGCGCTGGCGTCCAGTTC", "TGTGCTGGCGTCCAGTTC", "TGCGCAGGGGTCCAGTTC", "TGTGCAGGGGTCCAGTTC", "TGCGCCGGGGTCCAGTTC", "TGTGCCGGGGTCCAGTTC", "TGCGCGGGGGTCCAGTTC", "TGTGCGGGGGTCCAGTTC", "TGCGCTGGGGTCCAGTTC", "TGTGCTGGGGTCCAGTTC", "TGCGCAGGTGTCCAGTTC", "TGTGCAGGTGTCCAGTTC", "TGCGCCGGTGTCCAGTTC", "TGTGCCGGTGTCCAGTTC", "TGCGCGGGTGTCCAGTTC", "TGTGCGGGTGTCCAGTTC", "TGCGCTGGTGTCCAGTTC", "TGTGCTGGTGTCCAGTTC", "TGCGCAGGAGTGCAGTTC", "TGTGCAGGAGTGCAGTTC", "TGCGCCGGAGTGCAGTTC", "TGTGCCGGAGTGCAGTTC", "TGCGCGGGAGTGCAGTTC", "TGTGCGGGAGTGCAGTTC", "TGCGCTGGAGTGCAGTTC", "TGTGCTGGAGTGCAGTTC", "TGCGCAGGCGTGCAGTTC", "TGTGCAGGCGTGCAGTTC", "TGCGCCGGCGTGCAGTTC", "TGTGCCGGCGTGCAGTTC", "TGCGCGGGCGTGCAGTTC", "TGTGCGGGCGTGCAGTTC", "TGCGCTGGCGTGCAGTTC", "TGTGCTGGCGTGCAGTTC", "TGCGCAGGGGTGCAGTTC", "TGTGCAGGGGTGCAGTTC", "TGCGCCGGGGTGCAGTTC", "TGTGCCGGGGTGCAGTTC", "TGCGCGGGGGTGCAGTTC", "TGTGCGGGGGTGCAGTTC", "TGCGCTGGGGTGCAGTTC", "TGTGCTGGGGTGCAGTTC", "TGCGCAGGTGTGCAGTTC", "TGTGCAGGTGTGCAGTTC", "TGCGCCGGTGTGCAGTTC", "TGTGCCGGTGTGCAGTTC", "TGCGCGGGTGTGCAGTTC", "TGTGCGGGTGTGCAGTTC", "TGCGCTGGTGTGCAGTTC", "TGTGCTGGTGTGCAGTTC", "TGCGCAGGAGTTCAGTTC", "TGTGCAGGAGTTCAGTTC", "TGCGCCGGAGTTCAGTTC", "TGTGCCGGAGTTCAGTTC", "TGCGCGGGAGTTCAGTTC", "TGTGCGGGAGTTCAGTTC", "TGCGCTGGAGTTCAGTTC", "TGTGCTGGAGTTCAGTTC", "TGCGCAGGCGTTCAGTTC", "TGTGCAGGCGTTCAGTTC", "TGCGCCGGCGTTCAGTTC", "TGTGCCGGCGTTCAGTTC", "TGCGCGGGCGTTCAGTTC", "TGTGCGGGCGTTCAGTTC", "TGCGCTGGCGTTCAGTTC", "TGTGCTGGCGTTCAGTTC", "TGCGCAGGGGTTCAGTTC", "TGTGCAGGGGTTCAGTTC", "TGCGCCGGGGTTCAGTTC", "TGTGCCGGGGTTCAGTTC", "TGCGCGGGGGTTCAGTTC", "TGTGCGGGGGTTCAGTTC", "TGCGCTGGGGTTCAGTTC", "TGTGCTGGGGTTCAGTTC", "TGCGCAGGTGTTCAGTTC", "TGTGCAGGTGTTCAGTTC", "TGCGCCGGTGTTCAGTTC", "TGTGCCGGTGTTCAGTTC", "TGCGCGGGTGTTCAGTTC", "TGTGCGGGTGTTCAGTTC", "TGCGCTGGTGTTCAGTTC", "TGTGCTGGTGTTCAGTTC", "TGCGCAGGAGTACAATTT", "TGTGCAGGAGTACAATTT", "TGCGCCGGAGTACAATTT", "TGTGCCGGAGTACAATTT", "TGCGCGGGAGTACAATTT", "TGTGCGGGAGTACAATTT", "TGCGCTGGAGTACAATTT", "TGTGCTGGAGTACAATTT", "TGCGCAGGCGTACAATTT", "TGTGCAGGCGTACAATTT", "TGCGCCGGCGTACAATTT", "TGTGCCGGCGTACAATTT", "TGCGCGGGCGTACAATTT", "TGTGCGGGCGTACAATTT", "TGCGCTGGCGTACAATTT", "TGTGCTGGCGTACAATTT", "TGCGCAGGGGTACAATTT", "TGTGCAGGGGTACAATTT", "TGCGCCGGGGTACAATTT", "TGTGCCGGGGTACAATTT", "TGCGCGGGGGTACAATTT", "TGTGCGGGGGTACAATTT", "TGCGCTGGGGTACAATTT", "TGTGCTGGGGTACAATTT", "TGCGCAGGTGTACAATTT", "TGTGCAGGTGTACAATTT", "TGCGCCGGTGTACAATTT", "TGTGCCGGTGTACAATTT", "TGCGCGGGTGTACAATTT", "TGTGCGGGTGTACAATTT", "TGCGCTGGTGTACAATTT", "TGTGCTGGTGTACAATTT", "TGCGCAGGAGTCCAATTT", "TGTGCAGGAGTCCAATTT", "TGCGCCGGAGTCCAATTT", "TGTGCCGGAGTCCAATTT", "TGCGCGGGAGTCCAATTT", "TGTGCGGGAGTCCAATTT", "TGCGCTGGAGTCCAATTT", "TGTGCTGGAGTCCAATTT", "TGCGCAGGCGTCCAATTT", "TGTGCAGGCGTCCAATTT", "TGCGCCGGCGTCCAATTT", "TGTGCCGGCGTCCAATTT", "TGCGCGGGCGTCCAATTT", "TGTGCGGGCGTCCAATTT", "TGCGCTGGCGTCCAATTT", "TGTGCTGGCGTCCAATTT", "TGCGCAGGGGTCCAATTT", "TGTGCAGGGGTCCAATTT", "TGCGCCGGGGTCCAATTT", "TGTGCCGGGGTCCAATTT", "TGCGCGGGGGTCCAATTT", "TGTGCGGGGGTCCAATTT", "TGCGCTGGGGTCCAATTT", "TGTGCTGGGGTCCAATTT", "TGCGCAGGTGTCCAATTT", "TGTGCAGGTGTCCAATTT", "TGCGCCGGTGTCCAATTT", "TGTGCCGGTGTCCAATTT", "TGCGCGGGTGTCCAATTT", "TGTGCGGGTGTCCAATTT", "TGCGCTGGTGTCCAATTT", "TGTGCTGGTGTCCAATTT", "TGCGCAGGAGTGCAATTT", "TGTGCAGGAGTGCAATTT", "TGCGCCGGAGTGCAATTT", "TGTGCCGGAGTGCAATTT", "TGCGCGGGAGTGCAATTT", "TGTGCGGGAGTGCAATTT", "TGCGCTGGAGTGCAATTT", "TGTGCTGGAGTGCAATTT", "TGCGCAGGCGTGCAATTT", "TGTGCAGGCGTGCAATTT", "TGCGCCGGCGTGCAATTT", "TGTGCCGGCGTGCAATTT", "TGCGCGGGCGTGCAATTT", "TGTGCGGGCGTGCAATTT", "TGCGCTGGCGTGCAATTT", "TGTGCTGGCGTGCAATTT", "TGCGCAGGGGTGCAATTT", "TGTGCAGGGGTGCAATTT", "TGCGCCGGGGTGCAATTT", "TGTGCCGGGGTGCAATTT", "TGCGCGGGGGTGCAATTT", "TGTGCGGGGGTGCAATTT", "TGCGCTGGGGTGCAATTT", "TGTGCTGGGGTGCAATTT", "TGCGCAGGTGTGCAATTT", "TGTGCAGGTGTGCAATTT", "TGCGCCGGTGTGCAATTT", "TGTGCCGGTGTGCAATTT", "TGCGCGGGTGTGCAATTT", "TGTGCGGGTGTGCAATTT", "TGCGCTGGTGTGCAATTT", "TGTGCTGGTGTGCAATTT", "TGCGCAGGAGTTCAATTT", "TGTGCAGGAGTTCAATTT", "TGCGCCGGAGTTCAATTT", "TGTGCCGGAGTTCAATTT", "TGCGCGGGAGTTCAATTT", "TGTGCGGGAGTTCAATTT", "TGCGCTGGAGTTCAATTT", "TGTGCTGGAGTTCAATTT", "TGCGCAGGCGTTCAATTT", "TGTGCAGGCGTTCAATTT", "TGCGCCGGCGTTCAATTT", "TGTGCCGGCGTTCAATTT", "TGCGCGGGCGTTCAATTT", "TGTGCGGGCGTTCAATTT", "TGCGCTGGCGTTCAATTT", "TGTGCTGGCGTTCAATTT", "TGCGCAGGGGTTCAATTT", "TGTGCAGGGGTTCAATTT", "TGCGCCGGGGTTCAATTT", "TGTGCCGGGGTTCAATTT", "TGCGCGGGGGTTCAATTT", "TGTGCGGGGGTTCAATTT", "TGCGCTGGGGTTCAATTT", "TGTGCTGGGGTTCAATTT", "TGCGCAGGTGTTCAATTT", "TGTGCAGGTGTTCAATTT", "TGCGCCGGTGTTCAATTT", "TGTGCCGGTGTTCAATTT", "TGCGCGGGTGTTCAATTT", "TGTGCGGGTGTTCAATTT", "TGCGCTGGTGTTCAATTT", "TGTGCTGGTGTTCAATTT", "TGCGCAGGAGTACAGTTT", "TGTGCAGGAGTACAGTTT", "TGCGCCGGAGTACAGTTT", "TGTGCCGGAGTACAGTTT", "TGCGCGGGAGTACAGTTT", "TGTGCGGGAGTACAGTTT", "TGCGCTGGAGTACAGTTT", "TGTGCTGGAGTACAGTTT", "TGCGCAGGCGTACAGTTT", "TGTGCAGGCGTACAGTTT", "TGCGCCGGCGTACAGTTT", "TGTGCCGGCGTACAGTTT", "TGCGCGGGCGTACAGTTT", "TGTGCGGGCGTACAGTTT", "TGCGCTGGCGTACAGTTT", "TGTGCTGGCGTACAGTTT", "TGCGCAGGGGTACAGTTT", "TGTGCAGGGGTACAGTTT", "TGCGCCGGGGTACAGTTT", "TGTGCCGGGGTACAGTTT", "TGCGCGGGGGTACAGTTT", "TGTGCGGGGGTACAGTTT", "TGCGCTGGGGTACAGTTT", "TGTGCTGGGGTACAGTTT", "TGCGCAGGTGTACAGTTT", "TGTGCAGGTGTACAGTTT", "TGCGCCGGTGTACAGTTT", "TGTGCCGGTGTACAGTTT", "TGCGCGGGTGTACAGTTT", "TGTGCGGGTGTACAGTTT", "TGCGCTGGTGTACAGTTT", "TGTGCTGGTGTACAGTTT", "TGCGCAGGAGTCCAGTTT", "TGTGCAGGAGTCCAGTTT", "TGCGCCGGAGTCCAGTTT", "TGTGCCGGAGTCCAGTTT", "TGCGCGGGAGTCCAGTTT", "TGTGCGGGAGTCCAGTTT", "TGCGCTGGAGTCCAGTTT", "TGTGCTGGAGTCCAGTTT", "TGCGCAGGCGTCCAGTTT", "TGTGCAGGCGTCCAGTTT", "TGCGCCGGCGTCCAGTTT", "TGTGCCGGCGTCCAGTTT", "TGCGCGGGCGTCCAGTTT", "TGTGCGGGCGTCCAGTTT", "TGCGCTGGCGTCCAGTTT", "TGTGCTGGCGTCCAGTTT", "TGCGCAGGGGTCCAGTTT", "TGTGCAGGGGTCCAGTTT", "TGCGCCGGGGTCCAGTTT", "TGTGCCGGGGTCCAGTTT", "TGCGCGGGGGTCCAGTTT", "TGTGCGGGGGTCCAGTTT", "TGCGCTGGGGTCCAGTTT", "TGTGCTGGGGTCCAGTTT", "TGCGCAGGTGTCCAGTTT", "TGTGCAGGTGTCCAGTTT", "TGCGCCGGTGTCCAGTTT", "TGTGCCGGTGTCCAGTTT", "TGCGCGGGTGTCCAGTTT", "TGTGCGGGTGTCCAGTTT", "TGCGCTGGTGTCCAGTTT", "TGTGCTGGTGTCCAGTTT", "TGCGCAGGAGTGCAGTTT", "TGTGCAGGAGTGCAGTTT", "TGCGCCGGAGTGCAGTTT", "TGTGCCGGAGTGCAGTTT", "TGCGCGGGAGTGCAGTTT", "TGTGCGGGAGTGCAGTTT", "TGCGCTGGAGTGCAGTTT", "TGTGCTGGAGTGCAGTTT", "TGCGCAGGCGTGCAGTTT", "TGTGCAGGCGTGCAGTTT", "TGCGCCGGCGTGCAGTTT", "TGTGCCGGCGTGCAGTTT", "TGCGCGGGCGTGCAGTTT", "TGTGCGGGCGTGCAGTTT", "TGCGCTGGCGTGCAGTTT", "TGTGCTGGCGTGCAGTTT", "TGCGCAGGGGTGCAGTTT", "TGTGCAGGGGTGCAGTTT", "TGCGCCGGGGTGCAGTTT", "TGTGCCGGGGTGCAGTTT", "TGCGCGGGGGTGCAGTTT", "TGTGCGGGGGTGCAGTTT", "TGCGCTGGGGTGCAGTTT", "TGTGCTGGGGTGCAGTTT", "TGCGCAGGTGTGCAGTTT", "TGTGCAGGTGTGCAGTTT", "TGCGCCGGTGTGCAGTTT", "TGTGCCGGTGTGCAGTTT", "TGCGCGGGTGTGCAGTTT", "TGTGCGGGTGTGCAGTTT", "TGCGCTGGTGTGCAGTTT", "TGTGCTGGTGTGCAGTTT", "TGCGCAGGAGTTCAGTTT", "TGTGCAGGAGTTCAGTTT", "TGCGCCGGAGTTCAGTTT", "TGTGCCGGAGTTCAGTTT", "TGCGCGGGAGTTCAGTTT", "TGTGCGGGAGTTCAGTTT", "TGCGCTGGAGTTCAGTTT", "TGTGCTGGAGTTCAGTTT", "TGCGCAGGCGTTCAGTTT", "TGTGCAGGCGTTCAGTTT", "TGCGCCGGCGTTCAGTTT", "TGTGCCGGCGTTCAGTTT", "TGCGCGGGCGTTCAGTTT", "TGTGCGGGCGTTCAGTTT", "TGCGCTGGCGTTCAGTTT", "TGTGCTGGCGTTCAGTTT", "TGCGCAGGGGTTCAGTTT", "TGTGCAGGGGTTCAGTTT", "TGCGCCGGGGTTCAGTTT", "TGTGCCGGGGTTCAGTTT", "TGCGCGGGGGTTCAGTTT", "TGTGCGGGGGTTCAGTTT", "TGCGCTGGGGTTCAGTTT", "TGTGCTGGGGTTCAGTTT", "TGCGCAGGTGTTCAGTTT", "TGTGCAGGTGTTCAGTTT", "TGCGCCGGTGTTCAGTTT", "TGTGCCGGTGTTCAGTTT", "TGCGCGGGTGTTCAGTTT", "TGTGCGGGTGTTCAGTTT", "TGCGCTGGTGTTCAGTTT", "TGTGCTGGTGTTCAGTTT"};

    ClonotypeAA clonotype = aligner.buildClonotype();

    YMIR_ASSERT2(clonotype.nVar(), 3)
    YMIR_ASSERT2(clonotype.nDiv(), 2)
    YMIR_ASSERT2(clonotype.nJoi(), 3)

    MAAGaa maag = maag_builder.build(clonotype);

    std::vector<ClonotypeNuc> clonotype_vec;
    NaiveCDR3NucleotideAligner naligner(genes, VDJAlignerParameters(3,
                                                                    VDJAlignmentEventScore(AlignmentEventScore(1, -1, 1),
                                                                                           AlignmentEventScore(1, -1, 1),
                                                                                           AlignmentEventScore(1, -1, 1)),
                                                                    VDJAlignmentScoreThreshold(2, 2, 2)));
    for (auto &seq: rev_nuc) {
        naligner.setSequence(seq).setRecombination(VDJ_RECOMB);
        if (naligner.alignVar() && naligner.alignDiv() && naligner.alignJoi()) {
            clonotype_vec.push_back(naligner.buildClonotype());
        }
        naligner.clear();
    }

    prob_t sum_prob = 0;
    std::vector<prob_t> prob_vec;
    for (int i = 0; i < clonotype_vec.size(); ++i) {
        sum_prob += maag_builder.buildAndCompute(clonotype_vec[i], NO_ERRORS);
        auto tmp_maag = maag_builder.build(clonotype_vec[i], NO_METADATA, NO_ERRORS);
        for (int v_i = 0; v_i < tmp_maag.nVar(); ++v_i) {
            for (int d_i = 0; d_i < tmp_maag.nDiv(); ++d_i) {
                for (int j_i = 0; j_i < tmp_maag.nJoi(); ++j_i) {
                    prob_vec.push_back(tmp_maag.fullProbability(v_i, d_i, j_i));
                }
            }
        }
    }

    // TODO: need: 1.67487e-07
    // result: 4.32902e-08
    YMIR_ASSERT3(sum_prob, std::accumulate(prob_vec.begin(), prob_vec.end(), .0))
    YMIR_ASSERT3(maag.fullProbability(), std::accumulate(prob_vec.begin(), prob_vec.end(), .0))
    YMIR_ASSERT3(maag.fullProbability(), sum_prob)
    YMIR_ASSERT3(maag.fullProbability(), maag_builder.buildAndCompute(clonotype))

YMIR_TEST_END


YMIR_TEST_START(test_maag_forward_backward_vj)

    ModelParameterVector mvec = make_test_events_vj();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
//    alvec1.push_back("Vseg4");
    seqvec1.push_back("CCCA");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCCAGG");
//    seqvec1.push_back("TTCCCAGG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
//    alvec2.push_back("Jseg4");
    seqvec2.push_back("CCGTTT");
    seqvec2.push_back("CATT");
    seqvec2.push_back("AGGTTT");
//    seqvec2.push_back("AGGTTTGGG");

    VDJRecombinationGenes genes("VA", alvec1, seqvec1, "JA", alvec2, seqvec2);

    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)] = 1;
    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 1)] = 0;
    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 2)] = 0;
    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 3)] = 0;
    MAAGBuilder maag_builder(mvec, genes);

    ClonotypeNucBuilder cl_builder;
    // CCCAAAAAAATT
    //       CCGTTT
    //         CATT
    //       AGGTTT
    // CCCG.AC.GGTTT
    // poses:
    // vs: 0-1-2-3-4-5
    // js: 7-8-9-10-11
    cl_builder.setSequence("CCCAAAAAAATT")
            .addVarAlignment(1, 1, 1, 4)
            .addVarAlignment(3, 1, 1, 4)
//            .addVarAlignment(4, 3, 1, 4)
            .addJoiAlignment(1, 2, 11, 2)
            .addJoiAlignment(2, 2, 10, 3)
            .addJoiAlignment(3, 2, 11, 2);
//            .addJoiAlignment(4, 2, 11, 3);
    cl_builder.setRecombination(VJ_RECOMB);
    ClonotypeNuc clonotype = cl_builder.buildClonotype();

    MAAGnuc maag = maag_builder.build(clonotype, SAVE_METADATA, NO_ERRORS);

    MAAGForwardBackwardAlgorithm algo(maag);

    YMIR_ASSERT2(algo.status(), true)

//    while (!algo.is_empty()) {
//        auto temp = algo.nextEvent();
//        cout << temp.first << " : " << temp.second << endl;
//    }
//
//    cout << algo.VJ_nuc_probs()[0] << endl;
//    cout << algo.VJ_nuc_probs()[1] << endl;
//    cout << algo.VJ_nuc_probs()[2] << endl;
//    cout << algo.VJ_nuc_probs()[3] << endl;

    YMIR_ASSERT(abs(algo.fullProbability() - maag.fullProbability()) < 8e-20)

    YMIR_ASSERT(abs(algo.bfullProbability() - maag.fullProbability()) < 8e-20)

YMIR_TEST_END


YMIR_TEST_START(test_maag_forward_backward_vdj)

    ModelParameterVector mvec = make_test_events_vdj();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("CCCG");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCCGGG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
    seqvec2.push_back("CCGTTT");
    seqvec2.push_back("ATTT");
    seqvec2.push_back("AGGTTT");

    vector<string> alvec3;
    vector<string> seqvec3;
    alvec3.push_back("Dseg1");
    alvec3.push_back("Dseg2");
    alvec3.push_back("Dseg3");
    seqvec3.push_back("GTTT");
    seqvec3.push_back("ACCGGT");
    seqvec3.push_back("CCCGGAC");

    VDJRecombinationGenes genes("VB", alvec1, seqvec1, "JB", alvec2, seqvec2, "DB", alvec3, seqvec3);

    MAAGBuilder maag_builder(mvec, genes);

    ClonotypeNucBuilder cl_builder;
    /*
     D1:
       CCCGACGGTTT
             .GTTT
     D2:
       CCCGACGGTTT
      A.CCG.GT
         AC.CGGT

     D3:
       CCCGACGGTTT
     CCCG.GAC
       CCCG.GAC
         CC.CGG.AC
    */
    cl_builder.setSequence("CCCGACGGTTT")
            .addVarAlignment(1, 1, 1, 4)
            .addVarAlignment(3, 1, 1, 5)
            .addJoiAlignment(1, 2, 8, 4)
            .addJoiAlignment(2, 2, 9, 3)
            .addJoiAlignment(3, 2, 7, 5)
            .addDivAlignment(2, 2, 2, 3)
            .addDivAlignment(2, 3, 6, 4)
            .addDivAlignment(3, 5, 4, 3)
            .addDivAlignment(3, 1, 1, 4)
            .addDivAlignment(3, 3, 6, 3)
            .addDivAlignment(1, 1, 8, 4);
    cl_builder.setRecombination(VDJ_RECOMB);
    ClonotypeNuc clonotype = cl_builder.buildClonotype();

    MAAGnuc maag = maag_builder.build(clonotype, SAVE_METADATA, NO_ERRORS);

//    for (int node_i = 0; node_i < maag.chainSize(); ++node_i) {
//        for (int mat_i = 0; mat_i < maag.nodeSize(node_i); ++mat_i) {
//            for (int row_i = 0; row_i < maag.nodeRows(node_i); ++row_i) {
//                for (int col_i = 0; col_i < maag.nodeColumns(node_i); ++col_i) {
//                    cout << maag.event_index(node_i, mat_i, row_i, col_i) << endl;
//                }
//            }
//        }
//    }

    YMIR_ASSERT2(maag.nVar(), 2)
    YMIR_ASSERT2(maag.nDiv(), 3)
    YMIR_ASSERT2(maag.nJoi(), 3)

    MAAGForwardBackwardAlgorithm algo(maag);

    YMIR_ASSERT2(algo.status(), true)

    YMIR_ASSERT(abs(algo.fullProbability() - maag.fullProbability()) < 6e-20)

    YMIR_ASSERT(abs(algo.bfullProbability() - maag.fullProbability()) < 6e-20)

YMIR_TEST_END


int main(int argc, char* argv[]) {

    TEST_DATA_FOLDER = string(argv[1]) + string("/");

//    mpreal::set_default_prec(200);

    //
    // MEGA Todo: make good tests with some unit-testing framework (CTest)
    //

    //**************  INITIALISATION  **************//
    size_t tests_passed = 0, all_tests = 0;
    vector<TestInfo> failed_test_info;
    //**************  **************//



    //**************  TEST CASES  **************//

    // Test for Multi-Matrix Chains
//    YMIR_TEST(test_mmc())

    // Tests for MAAG / MAAG builder
    YMIR_TEST(test_maag_vj())
//    YMIR_TEST(test_maag_vj_err())
//    YMIR_TEST(test_maag_vdj())
//    YMIR_TEST(test_maag_vdj_err())
//    YMIR_TEST(test_maag_builder_replace_vj())
//    YMIR_TEST(test_maag_builder_replace_vdj())

//    YMIR_TEST(test_maag_vj_aa())
//    YMIR_TEST(test_maag_vdj_aa_very_simple())
//    YMIR_TEST(test_maag_vdj_aa_simple_D())
//    YMIR_TEST(test_maag_vdj_aa())

    // Tests for forward-backward algorithms
//    YMIR_TEST(test_maag_forward_backward_vj())
//    YMIR_TEST(test_maag_forward_backward_vdj())
    

    // Test for computing full nucleotide probabilities of repertoire with PAM.

    // Test for computing full amino acid probabilities of repertoire with PAM.

    // Tests for statistical inference of PAM parameters.

    //**************  **************//



    //**************  TESTING RESULTS  **************//
    std::cout << std::endl;
    cout << "Tests passed:\t" << tests_passed << endl;
    cout << "Tests failed:\t" << (all_tests - tests_passed) << endl;

    if (all_tests - tests_passed) {
        cout << "Failed tests:" << endl;
        for (size_t i = 0; i < failed_test_info.size(); ++i) {
            TestInfo ti = failed_test_info[i];
            cout << (i+1) << ":  " << ti.test_name << endl;
            for (size_t j = 0; j < ti.failed_cases.size(); ++j) {
                cout << '\t' << (int) (i+1) << '.' << (int) (j+1) << ":  " << ti.failed_cases[j] << endl;
            }
        }
    }
    //**************  **************//

    return 0;
}