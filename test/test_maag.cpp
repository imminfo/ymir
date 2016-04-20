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
    seqvec3.push_back("ACCGGT");
    seqvec3.push_back("CCCGGAC");

    VDJRecombinationGenes genes("VB", alvec1, seqvec1, "JB", alvec2, seqvec2, "DB", alvec3, seqvec3);

    MAAGBuilder maag_builder(mvec, genes);

    ClonotypeNucBuilder cl_builder;
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

//    cout << "V 0:" << maag.rows(0) << ":" << maag.cols(0) << endl;
//    cout << "Vdel 1:" << maag.rows(1) << ":" << maag.cols(1) << endl;
//    cout << "VDins 2:" << maag.rows(2) << ":" << maag.cols(2) << endl;
//    cout << "Ddel 3:" << maag.rows(3) << ":" << maag.cols(3) << endl;
//    cout << "DJins 4:" << maag.rows(4) << ":" << maag.cols(4) << endl;
//    cout << "Jdel 5:" << maag.rows(5) << ":" << maag.cols(5) << endl;
//    cout << "JD 6:" << maag.rows(6) << ":" << maag.cols(6) << endl;

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
    YMIR_ASSERT2(maag.event_index(4, 0, 0, 0), mvec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 3))
    YMIR_ASSERT2(maag.event_index(4, 0, 1, 1), mvec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 3))
    YMIR_ASSERT2(maag.event_index(4, 0, 3, 0), 0)

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

    ModelParameterVector mvec = make_test_events_vj4();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1"); seqvec1.push_back("TGTGC");
//    alvec1.push_back("Vseg2"); seqvec1.push_back("TGCG");
//    alvec1.push_back("Vseg3"); seqvec1.push_back("TGA");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1"); seqvec2.push_back("TTT");
//    alvec2.push_back("Jseg2"); seqvec2.push_back("ATTC");
//    alvec2.push_back("Jseg3"); seqvec2.push_back("AATT");

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
    CDR3AminoAcidAligner aligner(genes, VDJAlignerParameters(3));
    aligner.setSequence("C").setRecombination(VJ_RECOMB);
//    aligner.setSequence("CF").setRecombination(VJ_RECOMB);
//    aligner.setSequence("CAF").setRecombination(VJ_RECOMB);
//    aligner.setSequence("CASF").setRecombination(VJ_RECOMB);
    YMIR_ASSERT(aligner.alignVar())
    YMIR_ASSERT(aligner.alignJoi())

    ClonotypeAA clonotype = aligner.buildClonotype();

    YMIR_ASSERT2(clonotype.getVarLen(0), 5)
//    YMIR_ASSERT2(clonotype.getVarLen(1), 4)
//    YMIR_ASSERT2(clonotype.getVarLen(2), 2)
//    YMIR_ASSERT2(clonotype.getVarCodon(0, 1), compute_codon_hash({true, true, false, false, false, false}, 0))
//    YMIR_ASSERT2(clonotype.getVarCodon(0, 2), compute_codon_hash({true, true, false, false, false, false}, 0))
//    YMIR_ASSERT2(clonotype.getVarCodon(0, 3), compute_codon_hash({true, false, false, false, false, false}, 0))
//    YMIR_ASSERT2(clonotype.getVarCodon(0, 4), compute_codon_hash({true, true, true, true, false, false}, 0))
//    YMIR_ASSERT2(clonotype.getVarCodon(0, 5), compute_codon_hash({true, true, true, true, false, false}, 0))

    YMIR_ASSERT2(clonotype.getJoiLen(0), 3)
//    YMIR_ASSERT2(clonotype.getJoiLen(1), 4)
//    YMIR_ASSERT2(clonotype.getJoiLen(2), 2)
//    YMIR_ASSERT2(clonotype.getJoiSeqStart(2), 11)
//    YMIR_ASSERT2(clonotype.getJoiSeqEnd(2), 12)
//    YMIR_ASSERT2(clonotype.getJoiGeneStart(2), 3)
//    YMIR_ASSERT2(clonotype.getJoiGeneEnd(2), 4)
//    YMIR_ASSERT2(clonotype.getJoiCodon(0, 1), compute_codon_hash({true, false, false, false, false, false}, 0))
//    YMIR_ASSERT2(clonotype.getJoiCodon(0, 2), compute_codon_hash({true, false, false, false, false, false}, 0))
//    YMIR_ASSERT2(clonotype.getJoiCodon(0, 3), compute_codon_hash({true, false, false, false, false, false}, 0))

    MAAGaa maag = maag_builder.build(clonotype);

    YMIR_ASSERT2(maag.nVar(), 3)
    YMIR_ASSERT2(maag.nDiv(), 0)
    YMIR_ASSERT2(maag.nJoi(), 3)

    YMIR_ASSERT2(maag.n_poses(), 11)

//    YMIR_ASSERT2(maag.position(0), 0)
//    YMIR_ASSERT2(maag.position(1), 1)
//    YMIR_ASSERT2(maag.position(2), 2)
//    YMIR_ASSERT2(maag.position(3), 3)
//    YMIR_ASSERT2(maag.position(4), 4)
//    YMIR_ASSERT2(maag.position(5), 5)
//
//    YMIR_ASSERT2(maag.position(6), 9)
//    YMIR_ASSERT2(maag.position(7), 10)
//    YMIR_ASSERT2(maag.position(8), 11)
//    YMIR_ASSERT2(maag.position(9), 12)
//    YMIR_ASSERT2(maag.position(10), 13)

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

    rev_nuc = {"TGCGCATTC", "TGTGCATTC", "TGCGCCTTC", "TGTGCCTTC", "TGCGCGTTC", "TGTGCGTTC", "TGCGCTTTC", "TGTGCTTTC", "TGCGCATTT", "TGTGCATTT", "TGCGCCTTT", "TGTGCCTTT", "TGCGCGTTT", "TGTGCGTTT", "TGCGCTTTT", "TGTGCTTTT"};

    rev_nuc = { "TGCTTC", "TGTTTC", "TGCTTT", "TGTTTT" };

    rev_nuc = { "TGC", "TGT" };

    std::vector<ClonotypeNuc> clonotype_vec;
    NaiveCDR3NucleotideAligner naligner(genes, VDJAlignerParameters(3));
    for (auto &seq: rev_nuc) {
        naligner.setSequence(seq).setRecombination(VJ_RECOMB);
        YMIR_ASSERT(naligner.alignVar())
        YMIR_ASSERT(naligner.alignJoi())
        clonotype_vec.push_back(naligner.buildClonotype());
    }

    prob_t sum_prob = 0;
    for (int i = 0; i < clonotype_vec.size(); ++i) {
        sum_prob += maag_builder.buildAndCompute(clonotype_vec[i], NO_ERRORS);
//        std::cout << rev_nuc[i] << " = " << maag_builder.buildAndCompute(clonotype_vec[i], NO_ERRORS) << std::endl;
//        maag_builder.build(clonotype_vec[i], NO_METADATA, NO_ERRORS).print();
//        if (!maag_builder.buildAndCompute(clonotype_vec[i], NO_ERRORS)) {
//            maag_builder.build(clonotype_vec[i], NO_METADATA, NO_ERRORS).print();
//        }
    }

//    YMIR_ASSERT3(maag.fullProbability(0, 0), 0)
//    YMIR_ASSERT3(maag.fullProbability(1, 0), 0)
//    YMIR_ASSERT3(maag.fullProbability(2, 0), 0)
//    YMIR_ASSERT3(maag.fullProbability(0, 1), 0)
//    YMIR_ASSERT3(maag.fullProbability(1, 1), 0)
//    YMIR_ASSERT3(maag.fullProbability(2, 1), 0)
//    YMIR_ASSERT3(maag.fullProbability(0, 2), 0)
//    YMIR_ASSERT3(maag.fullProbability(1, 2), 0)
//    YMIR_ASSERT3(maag.fullProbability(2, 2), 0)

    YMIR_ASSERT3(maag.fullProbability(0, 0), 0)
    maag.print();
    maag.printCodons();

    YMIR_ASSERT3(maag.fullProbability(), sum_prob)
    YMIR_ASSERT3(maag.fullProbability(), maag_builder.buildAndCompute(clonotype))

YMIR_TEST_END


YMIR_TEST_START(test_maag_vdj_aa)

    YMIR_ASSERT(false)

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
    // MEGA TO-DO: make good tests with some unit-testing framework (CTest)
    //

    //**************  INITIALISATION  **************//
    size_t tests_passed = 0, all_tests = 0;
    vector<TestInfo> failed_test_info;
    //**************  **************//



    //**************  TEST CASES  **************//

    // Test for Multi-Matrix Chains
    YMIR_TEST(test_mmc())

    // Tests for MAAG / MAAG builder
    YMIR_TEST(test_maag_vj())
    YMIR_TEST(test_maag_vj_err())
    YMIR_TEST(test_maag_vdj())
    YMIR_TEST(test_maag_vdj_err())
    YMIR_TEST(test_maag_builder_replace_vj())
    YMIR_TEST(test_maag_builder_replace_vdj())

    YMIR_TEST(test_maag_vj_aa())
    YMIR_TEST(test_maag_vdj_aa())

    // Tests for forward-backward algorithms
    YMIR_TEST(test_maag_forward_backward_vj())
    YMIR_TEST(test_maag_forward_backward_vdj())
    

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