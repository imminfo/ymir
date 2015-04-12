/*
 * Ymir <imminfo.github.io/Ymir>
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


#define YMIR_TEST(res, s) {all_tests += 1; if (res.size() == 0) {tests_passed += 1;} \
                            else { \
                            failed_test_info.push_back(TestInfo(s, res)); \
                            } \
                            }

#define YMIR_TEST_START(funname) vector<string> funname() { vector<string> _failed_cases;
#define YMIR_ASSERT(expr) { if (!(expr)) { _failed_cases.push_back(#expr); } };
#define YMIR_ASSERT2(expr1, expr2) { if ((expr1) != (expr2)) { std::stringstream ss; ss << #expr1 << " == " << #expr2 << "  (result: " << (expr1) << ", need: " << (expr2) << ")";_failed_cases.push_back(ss.str()); } };
#define YMIR_TEST_END return _failed_cases; }

#define TEST_DATA_FOLDER string("/Users/vdn/ymir/test/data/")


#include <iostream>
#include <list>
#include <sstream>

#include "parser.h"
#include "statisticalinferencealgorithm.h"


using namespace std;
using namespace ymir;


YMIR_TEST_START(test_basic)
YMIR_ASSERT(1 == 1)
YMIR_TEST_END


YMIR_TEST_START(test_model_param_vec_vj)
    vector<prob_t> v1;  // param vec
    vector<eventind_t> v2;  // lens vec
    vector<eventind_t> v3;  // event classes
    vector<seq_len_t> v4;  // event family col numbers

    // V-J
    v1.push_back(.05); v1.push_back(.025); v1.push_back(.035); // J1
    v1.push_back(.045); v1.push_back(.055); v1.push_back(.065); // J2
    v1.push_back(.075); v1.push_back(.085); v1.push_back(.565); // J3
    v2.push_back(9);
    v3.push_back(0);
    v4.push_back(3);

    // V del
    v1.push_back(.75); v1.push_back(.25);
    v2.push_back(2);
    v4.push_back(0);

    v1.push_back(.4); v1.push_back(.5); v1.push_back(.1);
    v2.push_back(3);
    v4.push_back(0);

    v1.push_back(.3); v1.push_back(.1); v1.push_back(.2); v1.push_back(.4);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(1);

    // J del
    v1.push_back(.4); v1.push_back(.6);
    v2.push_back(2);
    v4.push_back(0);

    v1.push_back(.125); v1.push_back(.175); v1.push_back(.7);
    v2.push_back(3);
    v4.push_back(0);

    v3.push_back(4);
    v4.push_back(0);

    // VJ ins len
    v1.push_back(.31); v1.push_back(.39); v1.push_back(.1); v1.push_back(.1); v1.push_back(.1);
    v2.push_back(5);

    v3.push_back(6);
    v4.push_back(0);

    // VJ ins nuc
    // prev A
    v1.push_back(.05); v1.push_back(.08); v1.push_back(.03); v1.push_back(.84);
    v2.push_back(4);

    v3.push_back(7);
    v4.push_back(0);

    // prev C
    v1.push_back(.4); v1.push_back(.1); v1.push_back(.3); v1.push_back(.2);
    v2.push_back(4);

    v3.push_back(8);
    v4.push_back(0);

    // prev G
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2);
    v2.push_back(4);

    v3.push_back(9);
    v4.push_back(0);

    // prev T
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.25); v1.push_back(.3);
    v2.push_back(4);

    v3.push_back(10);
    v4.push_back(0);

    ModelParameterVector mvec(VJ_RECOMB, v1, v2, v3, v4);

    YMIR_ASSERT(mvec[0] == 0)

//    YMIR_ASSERT(mvec.prob_VJ_genes(1, 2) == .025)
//    YMIR_ASSERT(mvec.prob_VJ_genes(3, 1) == .075)
//
//    YMIR_ASSERT(mvec.prob_V_del(1, 0) == .75)
//    YMIR_ASSERT(mvec.prob_V_del(1, 1) == .25)
//    YMIR_ASSERT(mvec.prob_V_del(2, 0) == .4)
//    YMIR_ASSERT(mvec.prob_V_del(2, 1) == .5)
//    YMIR_ASSERT(mvec.prob_V_del(2, 2) == .1)
//    YMIR_ASSERT(mvec.prob_V_del(3, 0) == .3)
//    YMIR_ASSERT(mvec.prob_V_del(3, 3) == .4)
//
//    YMIR_ASSERT(mvec.prob_J_del(1, 0) == .4)
//    YMIR_ASSERT(mvec.prob_J_del(1, 1) == .6)
//    YMIR_ASSERT(mvec.prob_J_del(2, 0) == .125)
//    YMIR_ASSERT(mvec.prob_J_del(2, 2) == .7)
//
//    YMIR_ASSERT(mvec.prob_VJ_ins_len(0) - .31 < 1e-7)
//    YMIR_ASSERT(mvec.prob_VJ_ins_len(1) - .39 < 1e-7)
//    YMIR_ASSERT(mvec.prob_VJ_ins_len(2) - .1 < 1e-7)
//    YMIR_ASSERT(mvec.max_VJ_ins_len() == 4)

    YMIR_ASSERT2(mvec.eventClassSize(VJ_VAR_JOI_GEN), 9)
    YMIR_ASSERT2(mvec.eventClassSize(VJ_VAR_DEL), 9)
    YMIR_ASSERT2(mvec.eventFamilySize(VJ_VAR_DEL, 0), 2)
    YMIR_ASSERT2(mvec.eventFamilySize(VJ_VAR_DEL, 2), 4)
    YMIR_ASSERT2(mvec.eventClassSize(VJ_JOI_DEL), 5)
    YMIR_ASSERT2(mvec.eventFamilySize(VJ_JOI_DEL, 0), 2)
    YMIR_ASSERT2(mvec.eventFamilySize(VJ_JOI_DEL, 1), 3)
    YMIR_ASSERT2(mvec.eventClassSize(VJ_VAR_JOI_INS_LEN), 5)
    YMIR_ASSERT2(mvec.eventFamilySize(VJ_VAR_JOI_INS_LEN, 0), 5)
    YMIR_ASSERT2(mvec.eventClassSize(VJ_VAR_JOI_INS_NUC), 4)
    YMIR_ASSERT2(mvec.eventFamilySize(VJ_VAR_JOI_INS_NUC, 0), 4)

    YMIR_ASSERT2(mvec.event_prob(VJ_VAR_JOI_GEN, 0, 0, 1), .025)
    YMIR_ASSERT2(mvec.event_prob(VJ_VAR_JOI_GEN, 0, 2, 0), .075)

    YMIR_ASSERT2(mvec.event_prob(VJ_VAR_DEL, 0, 0), .75)
    YMIR_ASSERT2(mvec.event_prob(VJ_VAR_DEL, 0, 1), .25)
    YMIR_ASSERT2(mvec.event_prob(VJ_VAR_DEL, 1, 0), .4)
    YMIR_ASSERT2(mvec.event_prob(VJ_VAR_DEL, 1, 1), .5)
    YMIR_ASSERT2(mvec.event_prob(VJ_VAR_DEL, 1, 2), .1)
    YMIR_ASSERT2(mvec.event_prob(VJ_VAR_DEL, 2, 0), .3)
    YMIR_ASSERT2(mvec.event_prob(VJ_VAR_DEL, 2, 3), .4)

    YMIR_ASSERT2(mvec.event_prob(VJ_JOI_DEL, 0, 0), .4)
    YMIR_ASSERT2(mvec.event_prob(VJ_JOI_DEL, 0, 1), .6)
    YMIR_ASSERT2(mvec.event_prob(VJ_JOI_DEL, 1, 0), .125)
    YMIR_ASSERT2(mvec.event_prob(VJ_JOI_DEL, 1, 2), .7)

    YMIR_ASSERT2(round(mvec.event_prob(VJ_VAR_JOI_INS_LEN, 0, 0) * 100) / 100, .31)
    YMIR_ASSERT2(round(mvec.event_prob(VJ_VAR_JOI_INS_LEN, 0, 1) * 100) / 100, .39)
    YMIR_ASSERT2(round(mvec.event_prob(VJ_VAR_JOI_INS_LEN, 0, 2) * 100) / 100, .1)
    YMIR_ASSERT2(mvec.max_VJ_ins_len(), 4)

    YMIR_ASSERT2(mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)], .05)
    YMIR_ASSERT2(mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0) + 1], .08)
    YMIR_ASSERT2(mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0) + 4], .4)
YMIR_TEST_END


YMIR_TEST_START(test_model_param_vec_vdj)
    vector<prob_t> v1;  // param vec
    vector<eventind_t> v2;  // lens vec
    vector<eventind_t> v3;  // event classes
    vector<seq_len_t> v4;  // event family col numbers

    // V
    v1.push_back(.5); v1.push_back(.25); v1.push_back(.25);
    v2.push_back(3);
    v3.push_back(0);
    v4.push_back(0);

    // J-D (3 Js - 2 Ds)
    v1.push_back(.01); v1.push_back(.02);
    v1.push_back(.03); v1.push_back(.04);
    v1.push_back(.06); v1.push_back(.84);

    v2.push_back(6);

    v3.push_back(1);
    v4.push_back(2);

    // V del
    v1.push_back(.75); v1.push_back(.25);
    v2.push_back(2);
    v4.push_back(0);
    
    v1.push_back(.4); v1.push_back(.5); v1.push_back(.1);
    v2.push_back(3);
    v4.push_back(0);
    
    v1.push_back(.3); v1.push_back(.1); v1.push_back(.2); v1.push_back(.4);
    v2.push_back(4);
    v4.push_back(0);
    
    v3.push_back(2);    

    // J del
    v1.push_back(.4); v1.push_back(.6);
    v2.push_back(2);
    v4.push_back(0);

    v1.push_back(.125); v1.push_back(.175); v1.push_back(.7);
    v2.push_back(3);
    v4.push_back(0);

    v3.push_back(5);

    // D1 dels
    v1.push_back(.17); v1.push_back(.27);
    v1.push_back(.37); v1.push_back(.19);

    v2.push_back(4);
    v4.push_back(2);

    // D2 dels
    // 3 rows 2 columns
    v1.push_back(.11); v1.push_back(.12);
    v1.push_back(.13); v1.push_back(.14);
    v1.push_back(.15); v1.push_back(.35);

    v2.push_back(6);
    v4.push_back(2);

    v3.push_back(7);

    // VD ins len
    v1.push_back(.76); v1.push_back(.24);
    v2.push_back(2);
    v4.push_back(0);

    v3.push_back(9);

    // DJ ins len
    v1.push_back(.89); v1.push_back(.10); v1.push_back(.01);
    v2.push_back(3);
    v4.push_back(0);

    v3.push_back(10);

    // VD ins nuc
    // prev A
    v1.push_back(.05); v1.push_back(.08); v1.push_back(.03); v1.push_back(.84);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(11);

    // prev C
    v1.push_back(.4); v1.push_back(.1); v1.push_back(.3); v1.push_back(.2);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(12);

    // prev G
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(13);

    // prev T
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.25); v1.push_back(.3);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(14);

    // DJ ins nuc
    // prev A
    v1.push_back(.009); v1.push_back(.06); v1.push_back(.48); v1.push_back(.451);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(15);

    // prev C
    v1.push_back(.39); v1.push_back(.01); v1.push_back(.31); v1.push_back(.29);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(16);

    // prev G
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(17);

    // prev T
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.25); v1.push_back(.3);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(18);

    ModelParameterVector mvec(VDJ_RECOMB, v1, v2, v3, v4);

    YMIR_ASSERT2(mvec[0], 0)

//    YMIR_ASSERT(mvec.prob_V_gene(1) == .5)
//    YMIR_ASSERT(mvec.prob_V_gene(3) == .25)
//
//    YMIR_ASSERT(mvec.prob_JD_genes(1, 1) == .01)
//    YMIR_ASSERT(mvec.prob_JD_genes(1, 2) == .02)
//    YMIR_ASSERT(mvec.prob_JD_genes(3, 1) == .06)
//    YMIR_ASSERT(mvec.prob_JD_genes(3, 2) == .84)
//
//    YMIR_ASSERT(mvec.prob_V_del(1, 0) == .75)
//    YMIR_ASSERT(mvec.prob_V_del(1, 1) == .25)
//    YMIR_ASSERT(mvec.prob_V_del(2, 0) == .4)
//    YMIR_ASSERT(mvec.prob_V_del(2, 1) == .5)
//    YMIR_ASSERT(mvec.prob_V_del(2, 2) == .1)
//    YMIR_ASSERT(mvec.prob_V_del(3, 0) == .3)
//    YMIR_ASSERT(mvec.prob_V_del(3, 3) == .4)
//
//    YMIR_ASSERT(mvec.prob_J_del(1, 0) == .4)
//    YMIR_ASSERT(mvec.prob_J_del(1, 1) == .6)
//    YMIR_ASSERT(mvec.prob_J_del(2, 0) == .125)
//    YMIR_ASSERT(mvec.prob_J_del(2, 2) == .7)
//
//    YMIR_ASSERT(mvec.prob_D_del(1, 0, 0) == .17)
//    YMIR_ASSERT(mvec.prob_D_del(1, 0, 1) == .27)
//    YMIR_ASSERT(mvec.prob_D_del(1, 1, 0) == .37)
//    YMIR_ASSERT(mvec.prob_D_del(1, 1, 1) == .19)
//
//    YMIR_ASSERT(mvec.prob_D_del(2, 0, 0) == .11)
//    YMIR_ASSERT(mvec.prob_D_del(2, 0, 1) == .12)
//    YMIR_ASSERT(mvec.prob_D_del(2, 2, 0) == .15)
//    YMIR_ASSERT(mvec.prob_D_del(2, 2, 1) == .35)
//
//    YMIR_ASSERT(mvec.prob_VD_ins_len(0) == .76)
//    YMIR_ASSERT(mvec.prob_VD_ins_len(1) == .24)
//    YMIR_ASSERT(mvec.max_VD_ins_len() == 1)
//
//    YMIR_ASSERT(mvec.prob_DJ_ins_len(0) == .89)
//    YMIR_ASSERT(mvec.prob_DJ_ins_len(1) == .10)
//    YMIR_ASSERT(mvec.prob_DJ_ins_len(2) == .01)
//    YMIR_ASSERT(mvec.max_DJ_ins_len() == 2)
//
//    YMIR_ASSERT(mvec[mvec.index_VD_ins_nuc()] == .05)
//    YMIR_ASSERT(mvec[mvec.index_VD_ins_nuc() + 1] == .08)
//    YMIR_ASSERT(mvec[mvec.index_VD_ins_nuc() + 4] == .4)
//
//    YMIR_ASSERT(mvec[mvec.index_DJ_ins_nuc()] == .009)
//    YMIR_ASSERT(mvec[mvec.index_DJ_ins_nuc() + 1] == .06)
//    YMIR_ASSERT(mvec[mvec.index_DJ_ins_nuc() + 4] == .39)
//
//    YMIR_ASSERT(mvec[mvec.index_VJ_ins_nuc()] == .05)
//    YMIR_ASSERT(mvec[mvec.index_VJ_ins_nuc() + 1] == .08)
//    YMIR_ASSERT(mvec[mvec.index_VJ_ins_nuc() + 4] == .4)

    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_GEN, 0, 0), .5)
    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_GEN, 0, 2), .25)

    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 0, 0), .01)
    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 0, 1), .02)
    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 1, 0), .03)
    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 1, 1), .04)
    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 2, 0), .06)
    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 2, 1), .84)

    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_DEL, 0, 0), .75)
    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_DEL, 0, 1), .25)
    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_DEL, 1, 0), .4)
    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_DEL, 1, 1), .5)
    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_DEL, 1, 2), .1)
    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_DEL, 2, 0), .3)
    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_DEL, 2, 3), .4)

    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DEL, 0, 0), .4)
    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DEL, 0, 1), .6)
    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DEL, 1, 0), .125)
    YMIR_ASSERT2(mvec.event_prob(VDJ_JOI_DEL, 1, 2), .7)

    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_DEL, 0, 0, 0), .17)
    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_DEL, 0, 0, 1), .27)
    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_DEL, 0, 1, 0), .37)
    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_DEL, 0, 1, 1), .19)

    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_DEL, 1, 0, 0), .11)
    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_DEL, 1, 0, 1), .12)
    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_DEL, 1, 2, 0), .15)
    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_DEL, 1, 2, 1), .35)

    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_DIV_INS_LEN, 0, 0), .76)
    YMIR_ASSERT2(mvec.event_prob(VDJ_VAR_DIV_INS_LEN, 0, 1), .24)
    YMIR_ASSERT2(mvec.max_VD_ins_len(), 1)

    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_JOI_INS_LEN, 0, 0), .89)
    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_JOI_INS_LEN, 0, 1), .10)
    YMIR_ASSERT2(mvec.event_prob(VDJ_DIV_JOI_INS_LEN, 0, 2), .01)
    YMIR_ASSERT2(mvec.max_DJ_ins_len(), 2)

    YMIR_ASSERT2(mvec[mvec.event_index(VDJ_VAR_DIV_INS_NUC, 0, 0)], .05)
    YMIR_ASSERT2(mvec[mvec.event_index(VDJ_VAR_DIV_INS_NUC, 0, 0) + 1], .08)
    YMIR_ASSERT2(mvec[mvec.event_index(VDJ_VAR_DIV_INS_NUC, 0, 0) + 4], .4)

    YMIR_ASSERT2(mvec[mvec.event_index(VDJ_DIV_JOI_INS_NUC, 0, 0)], .009)
    YMIR_ASSERT2(mvec[mvec.event_index(VDJ_DIV_JOI_INS_NUC, 0, 0) + 1], .06)
    YMIR_ASSERT2(mvec[mvec.event_index(VDJ_DIV_JOI_INS_NUC, 0, 0) + 4], .39)
YMIR_TEST_END


YMIR_TEST_START(test_genesegmentalphabet)
    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("ACT");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCC");

    GeneSegmentAlphabet gsa("testseg", alvec1, seqvec1);

    YMIR_ASSERT(gsa.name() == "testseg")
    YMIR_ASSERT(gsa.size() == 4)
    YMIR_ASSERT(gsa[0].sequence.size() == gsa["other"].sequence.size())
    YMIR_ASSERT(gsa[1].sequence == "ACT")
    YMIR_ASSERT(gsa[3].sequence == gsa["Vseg3"].sequence)

    gsa.appendPalindromicNucleotides(2, 0);
    YMIR_ASSERT(gsa[1].sequence == "GTACT")

    gsa.appendPalindromicNucleotides(0, 2);
    YMIR_ASSERT(gsa[1].sequence == "GTACTAG")

    gsa.appendPalindromicNucleotides(2, 2);
    YMIR_ASSERT(gsa[1].sequence == "CAGTACTAGCT")
YMIR_TEST_END


YMIR_TEST_START(test_vdjgenes2)
    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("ACT");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCC");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
    seqvec2.push_back("GGG");
    seqvec2.push_back("TTT");
    seqvec2.push_back("AAA");

    VDJRecombinationGenes vdjgens("VA", alvec1, seqvec1, "JA", alvec2, seqvec2);

    YMIR_ASSERT(!vdjgens.is_vdj())
    YMIR_ASSERT(vdjgens.V().size() == 4)
    YMIR_ASSERT(vdjgens.J().size() == 4)
    YMIR_ASSERT2(vdjgens.V()["RANDOM"].index, 0)
    YMIR_ASSERT2(vdjgens.V()[100].index, 0)
    YMIR_ASSERT2(vdjgens.J()["RANDOM2"].index, 0)
    YMIR_ASSERT2(vdjgens.J()[100].index, 0)
    YMIR_ASSERT(vdjgens.V()["Vseg1"].sequence == vdjgens.V()[1].sequence)
    YMIR_ASSERT(vdjgens.J()["Jseg2"].sequence == vdjgens.J()[2].sequence)
YMIR_TEST_END


YMIR_TEST_START(test_vdjgenes3)
    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    seqvec1.push_back("ACT");
    seqvec1.push_back("GGG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    seqvec2.push_back("GGG");
    seqvec2.push_back("TTT");

    vector<string> alvec3;
    vector<string> seqvec3;
    alvec3.push_back("Dseg1");
    alvec3.push_back("Dseg2");
    alvec3.push_back("Dseg3");
    seqvec3.push_back("GGAC");
    seqvec3.push_back("TTTA");
    seqvec3.push_back("AAAC");

    VDJRecombinationGenes vdjgens("VA", alvec1, seqvec1, "JA", alvec2, seqvec2, "DA", alvec3, seqvec3);

    YMIR_ASSERT(vdjgens.is_vdj())
    YMIR_ASSERT(vdjgens.V().size() == 3)
    YMIR_ASSERT(vdjgens.J().size() == 3)
    YMIR_ASSERT(vdjgens.D().size() == 4)
    YMIR_ASSERT(vdjgens.V()["Vseg1"].sequence == vdjgens.V()[1].sequence)
    YMIR_ASSERT(vdjgens.J()["Jseg2"].sequence == vdjgens.J()[2].sequence)
    YMIR_ASSERT(vdjgens.D()["Dseg3"].sequence == "AAAC")
YMIR_TEST_END


YMIR_TEST_START(test_genesegmentalphabet_read)
    // assert read
    bool ok;
    GeneSegmentAlphabet gsa_n("testseg", TEST_DATA_FOLDER + "RANDOM_FILE.txt", &ok);
    YMIR_ASSERT(!ok)
    GeneSegmentAlphabet gsa("testseg", TEST_DATA_FOLDER + "vgene.txt", &ok);
    YMIR_ASSERT(ok)

    YMIR_ASSERT(gsa.name() == "testseg")
    YMIR_ASSERT(gsa.size() == 4)
    YMIR_ASSERT(gsa[0].sequence.size() == gsa["other"].sequence.size())
    YMIR_ASSERT(gsa[1].sequence == "ACT")
    YMIR_ASSERT(gsa[3].sequence == gsa["Vseg3"].sequence)


    // assert write and than read again
    YMIR_ASSERT(gsa.write(TEST_DATA_FOLDER + "vgene_towrite.txt"))

    GeneSegmentAlphabet gsa1("testseg", TEST_DATA_FOLDER + "vgene_towrite.txt");

    YMIR_ASSERT(gsa1.name() == "testseg")
    YMIR_ASSERT(gsa1.size() == 4)
    YMIR_ASSERT(gsa1[0].sequence.size() == gsa1["other"].sequence.size())
    YMIR_ASSERT(gsa1[1].sequence == "ACT")
    YMIR_ASSERT(gsa1[3].sequence == gsa1["Vseg3"].sequence)
YMIR_TEST_END


YMIR_TEST_START(test_vdjgenes_read)
    VDJRecombinationGenes vdjgens("Vgene", TEST_DATA_FOLDER + "vgene.txt"
            , "Jgene", TEST_DATA_FOLDER + "jgene.txt"
            , "Dgene", TEST_DATA_FOLDER + "dgene.txt");

    // assert read
    YMIR_ASSERT(vdjgens.is_vdj())
    YMIR_ASSERT(vdjgens.V().size() == 4)
    YMIR_ASSERT(vdjgens.J().size() == 3)
    YMIR_ASSERT(vdjgens.D().size() == 4)
    YMIR_ASSERT(vdjgens.V()["Vseg1"].sequence == vdjgens.V()[1].sequence)
    YMIR_ASSERT(vdjgens.J()["Jseg2"].sequence == vdjgens.J()[2].sequence)
    YMIR_ASSERT(vdjgens.D()["Dseg3"].sequence == "AAAC")

    // assert write and than read again
    YMIR_ASSERT(vdjgens.write(TEST_DATA_FOLDER + "vgene_towrite.txt"
            , TEST_DATA_FOLDER + "jgene_towrite.txt"
            , TEST_DATA_FOLDER + "dgene_towrite.txt"));

    VDJRecombinationGenes vdjgens1("Vgene", TEST_DATA_FOLDER + "vgene_towrite.txt"
            , "Jgene", TEST_DATA_FOLDER + "jgene_towrite.txt"
            , "Dgene", TEST_DATA_FOLDER + "dgene_towrite.txt");

    YMIR_ASSERT(vdjgens1.is_vdj())
    YMIR_ASSERT(vdjgens1.V().size() == 4)
    YMIR_ASSERT(vdjgens1.J().size() == 3)
    YMIR_ASSERT(vdjgens1.D().size() == 4)
    YMIR_ASSERT(vdjgens1.V()["Vseg1"].sequence == vdjgens1.V()[1].sequence)
    YMIR_ASSERT(vdjgens1.J()["Jseg2"].sequence == vdjgens1.J()[2].sequence)
    YMIR_ASSERT(vdjgens1.D()["Dseg3"].sequence == "AAAC")
YMIR_TEST_END


YMIR_TEST_START(test_clone)
    segindex_t *segs = new segindex_t[10];
    segs[0] = 3;
    segs[1] = 2;
    segs[2] = 2;

    segs[3] = 1;
    segs[4] = 2;
    segs[5] = 3;

    segs[6] = 4;
    segs[7] = 5;

    segs[8] = 6;
    segs[9] = 7;

    seq_len_t *alignments = new seq_len_t[3 + 2 + 2*4 + 2*4];
    // V
    alignments[0] = 15;
    alignments[1] = 20;
    alignments[2] = 12;
    // J
    alignments[3] = 8;
    alignments[4] = 7;
    // 2 X 4tuples for D1
    alignments[5] = 1;
    alignments[6] = 2;
    alignments[7] = 3;
    alignments[8] = 4;

    alignments[9] = 1;
    alignments[10] = 2;
    alignments[11] = 6;
    alignments[12] = 7;

    // 2 X 4tuples for D2
    alignments[13] = 11;
    alignments[14] = 12;
    alignments[15] = 13;
    alignments[16] = 14;

    alignments[17] = 20;
    alignments[18] = 21;
    alignments[19] = 22;
    alignments[20] = 23;

    seq_len_t *nd = new seq_len_t[2];
    nd[0] = 2;
    nd[1] = 2;

    Clonotype c("cloneseq", true, segs, alignments, nd);
    YMIR_ASSERT(c.is_vdj())
    YMIR_ASSERT(c.getVar(0) == 1)
    YMIR_ASSERT(c.getVar(1) == 2)
    YMIR_ASSERT(c.getVar(2) == 3)
    YMIR_ASSERT(c.getJoi(0) == 4)
    YMIR_ASSERT(c.getJoi(1) == 5)
    YMIR_ASSERT(c.getDiv(0) == 6)
    YMIR_ASSERT(c.getDiv(1) == 7)
    YMIR_ASSERT(c.nVar() == 3)
    YMIR_ASSERT(c.nJoi() == 2)
    YMIR_ASSERT(c.nDiv() == 2)
    YMIR_ASSERT(c.nDalignments(0) == 2)
    YMIR_ASSERT(c.getVend(0) == 15)
    YMIR_ASSERT(c.getVend(2) == 12)
    YMIR_ASSERT(c.getJstart(1) == 7)
    YMIR_ASSERT(c.getDalignment(0, 0).Dstart == 1)
    YMIR_ASSERT(c.getDalignment(0, 0).Dend == 2)
    YMIR_ASSERT(c.getDalignment(0, 0).seqstart == 3)
    YMIR_ASSERT(c.getDalignment(0, 0).seqend == 4)
    YMIR_ASSERT(c.getDalignment(0, 1).seqstart == 6)
    YMIR_ASSERT(c.getDalignment(0, 1).seqend == 7)
    YMIR_ASSERT(c.getDalignment(1, 0).Dstart == 11)
    YMIR_ASSERT(c.getDalignment(1, 0).Dend == 12)
    YMIR_ASSERT(c.getDalignment(1, 0).seqstart == 13)
    YMIR_ASSERT(c.getDalignment(1, 0).seqend == 14)
    YMIR_ASSERT(c.getDalignment(1, 1).Dstart == 20)
    YMIR_ASSERT(c.getDalignment(1, 1).Dend == 21)
    YMIR_ASSERT(c.getDalignment(1, 1).seqstart == 22)
    YMIR_ASSERT(c.getDalignment(1, 1).seqend == 23)
YMIR_TEST_END


YMIR_TEST_START(test_clonebuilder_clonealign)
    ClonotypeBuilder cb;

    cb.setNucleotideSeq();
    cb.setSequence("nuclseq");
    cb.addValignment(10, 15)
            .addValignment(11, 25)
            .addValignment(12, 35)
            .addJalignment(20, 21)
            .addDalignment(31, 8, 9, 11, 12)
            .addDalignment(31, 8, 9, 13, 14)
            .addDalignment(30, 1, 2, 3, 4)
            .addDalignment(30, 1, 2, 5, 6)
            .addDalignment(32, 1, 2, 5, 6);

    Clonotype c = cb.buildClonotype();

    YMIR_ASSERT(c.sequence() == "nuclseq")
    YMIR_ASSERT(c.is_nucleotide())
    YMIR_ASSERT(c.is_vdj())
    YMIR_ASSERT(c.getVar(0) == 10)
    YMIR_ASSERT(c.getVend(0) == 15)
    YMIR_ASSERT(c.getVar(1) == 11)
    YMIR_ASSERT(c.getVend(1) == 25)
    YMIR_ASSERT(c.getVar(2) == 12)
    YMIR_ASSERT(c.getVend(2) == 35)
    YMIR_ASSERT(c.getJoi(0) == 20)
    YMIR_ASSERT(c.getJstart(0) == 21)
    YMIR_ASSERT(c.getDiv(0) == 31)
    YMIR_ASSERT(c.getDalignment(0, 0) == d_alignment_t(8,9,11,12))
    YMIR_ASSERT(c.getDalignment(0, 1) == d_alignment_t(8,9,13,14))
    YMIR_ASSERT(c.getDiv(1) == 30)
    YMIR_ASSERT(c.getDalignment(1, 0) == d_alignment_t(1,2,3,4))
    YMIR_ASSERT(c.getDalignment(1, 1) == d_alignment_t(1,2,5,6))
    YMIR_ASSERT(c.getDiv(2) == 32)
    YMIR_ASSERT(c.getDalignment(2, 0) == d_alignment_t(1,2,5,6))
    YMIR_ASSERT(c.nVar() == 3)
    YMIR_ASSERT(c.nJoi() == 1)
    YMIR_ASSERT(c.nDiv() == 3)
YMIR_TEST_END


YMIR_TEST_START(test_nuc_aligner)
    NaiveNucleotideAligner nna;
    YMIR_ASSERT(nna.align5end("ACGT", "ACGTT") == 4)
    YMIR_ASSERT(nna.align5end("ACGT", "ACGT") == 4)
    YMIR_ASSERT(nna.align5end("ACGT", "ACG") == 3)
    YMIR_ASSERT(nna.align5end("ACGT", "TTT") == 0)

    YMIR_ASSERT(nna.align3end("ACGT", "CGT") == 3)
    YMIR_ASSERT(nna.align3end("ACGT", "TACGT") == 4)
    YMIR_ASSERT(nna.align3end("ACGT", "TTCGT") == 3)
    YMIR_ASSERT(nna.align3end("ACGG", "TTTTT") == 0)


    YMIR_ASSERT(nna.alignLocal("AA", "TTAATAA", 3).size() == 0)
    YMIR_ASSERT(nna.alignLocal("AA", "TTAATAA", 2).size() == 2)
    YMIR_ASSERT(nna.alignLocal("AACCTT", "AAGGTTGGGGGTT", 2).size() == 3)
    YMIR_ASSERT(nna.alignLocal("ACT", "ACTGACGACGGTATCTAC", 2).size() == 5)
YMIR_TEST_END


YMIR_TEST_START(test_aa_aligner)
    YMIR_ASSERT(false)
YMIR_TEST_END


YMIR_TEST_START(test_mitcr_vj)
    RepertoireParser parser;
    YMIR_ASSERT(parser.loadConfig(TEST_DATA_FOLDER + "../../parsers/mitcr.json"))

    bool V_err, J_err;
    VDJRecombinationGenes vdj_genes("Vgene", TEST_DATA_FOLDER + "vgene.real.txt"
            , "Jgene", TEST_DATA_FOLDER + "jgene.real.txt", &V_err, &J_err);
    YMIR_ASSERT(V_err)
    YMIR_ASSERT(J_err)

    Cloneset cr;
    YMIR_ASSERT(parser.parse(TEST_DATA_FOLDER + "mitcr.alpha.txt", &cr, vdj_genes))

    YMIR_ASSERT(cr.size() == 30)
    YMIR_ASSERT(cr[0].sequence() == "TGTGCAGCAAGTACCCCCTTAAGCTGGTGGTACTAGCTATGGAAAGCTGACATTT")
    YMIR_ASSERT(!cr[0].is_vdj())
    YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(0)].allele == "TRAV13-1")
    YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(1)].allele == "TRAV13-2")
    YMIR_ASSERT(cr[0].nVar() == 2)
    YMIR_ASSERT(vdj_genes.J()[cr[0].getJoi(0)].allele == "TRAJ52")
    YMIR_ASSERT(cr[0].nJoi() == 1)

    YMIR_ASSERT(cr[2].sequence() == "TGTGCAACTCTTAGCAGGGATGAACACAGGCTTTCAGAAACTTGTATTT")
    YMIR_ASSERT(!cr[2].is_vdj())
    YMIR_ASSERT(vdj_genes.V()[cr[2].getVar(0)].allele == "TRAV12-3")
    YMIR_ASSERT(cr[2].nVar() == 1)
    YMIR_ASSERT(vdj_genes.J()[cr[2].getJoi(0)].allele == "TRAJ8")
    YMIR_ASSERT(vdj_genes.J()[cr[2].getJoi(1)].allele == "TRAJ18")
    YMIR_ASSERT(cr[2].nJoi() == 2)
YMIR_TEST_END


YMIR_TEST_START(test_mitcr_vdj)
    YMIR_ASSERT(false)
YMIR_TEST_END


YMIR_TEST_START(test_clorep)
    /*
    [index]
    [vector]
    head
    slice
    sample w/ replacement
    sample w/o replacement
     */
    YMIR_ASSERT(false)
YMIR_TEST_END


YMIR_TEST_START(test_clorep_view)
    YMIR_ASSERT(false)
YMIR_TEST_END


YMIR_TEST_START(test_markovchain_nuc)
    event_matrix_t mat;
    mat.resize(4, 4);
    // A
    mat(0, 0) = .1;
    mat(1, 0) = .4;
    mat(2, 0) = .25;
    mat(3, 0) = .25;
    // C
    mat(0, 1) = .7;
    mat(1, 1) = .1;
    mat(2, 1) = .1;
    mat(3, 1) = .1;
    // G
    mat(0, 2) = .3;
    mat(1, 2) = .3;
    mat(2, 2) = .15;
    mat(3, 2) = .25;
    // T
    mat(0, 3) = .1;
    mat(1, 3) = .4;
    mat(2, 3) = .2;
    mat(3, 3) = .3;

    vector<prob_t> vec;
    // prev A
    vec.push_back(.1);
    vec.push_back(.7);
    vec.push_back(.3);
    vec.push_back(.1);
    // prev C
    vec.push_back(.4);
    vec.push_back(.1);
    vec.push_back(.3);
    vec.push_back(.4);
    // prev G
    vec.push_back(.25);
    vec.push_back(.1);
    vec.push_back(.15);
    vec.push_back(.2);
    // prev T
    vec.push_back(.25);
    vec.push_back(.1);
    vec.push_back(.25);
    vec.push_back(.3);

    MarkovChain mc(mat);
    string s = "ACGT";

    // .7 * .3 * .2 = .042
    YMIR_ASSERT(mc.nucProbability("") == 0);
    YMIR_ASSERT(mc.nucProbability(s) == .042);
    YMIR_ASSERT(mc.nucProbability(s.begin(), 4) == .042);
    YMIR_ASSERT(mc.nucProbability(s.begin() + 1, 3) == .06);
    YMIR_ASSERT(mc.nucProbability(s.begin(), 0) == 0);

    mc.updateProbabilities(vec.begin());

    YMIR_ASSERT(mc.nucProbability(s) == .042);
    YMIR_ASSERT(mc.nucProbability(s.begin(), 4) == .042);
YMIR_TEST_END


YMIR_TEST_START(test_markovchain_aa)
    YMIR_ASSERT(false)
YMIR_TEST_END


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
            * mat.matrix(4, 0))(0, 0) == 4.4)

    YMIR_ASSERT((mat.matrix(0, 1)
            * mat.matrix(1, 1)
            * mat.matrix(2, 0)
            * mat.matrix(3, 1)
            * mat.matrix(4, 1))(0, 0) - 52.8 < 1e-14)
YMIR_TEST_END


YMIR_TEST_START(test_maag_vj)
    vector<prob_t> v1;  // param vec
    vector<eventind_t> v2;  // lens vec
    vector<eventind_t> v3;  // event classes
    vector<seq_len_t> v4;  // event family col numbers

    // V-J
    v1.push_back(.05); v1.push_back(.025); v1.push_back(.035); // J1
    v1.push_back(.045); v1.push_back(.055); v1.push_back(.065); // J2
    v1.push_back(.075); v1.push_back(.085); v1.push_back(.565); // J3
    v2.push_back(9);

    v3.push_back(0);
    v4.push_back(3);

    // V del
    v1.push_back(.4); v1.push_back(.5); v1.push_back(.05); v1.push_back(.02); v1.push_back(.03);
    v2.push_back(5);

    v1.push_back(.3); v1.push_back(.1); v1.push_back(.2); v1.push_back(.4);
    v2.push_back(4);

    v1.push_back(.75); v1.push_back(.005); v1.push_back(.01); v1.push_back(.02);
    v1.push_back(.03); v1.push_back(.04); v1.push_back(.145);
    v2.push_back(7);

    v3.push_back(1);
    v4.push_back(0);
    v4.push_back(0);
    v4.push_back(0);

    // J del
    v1.push_back(.07); v1.push_back(.2); v1.push_back(.3);
    v1.push_back(.34); v1.push_back(.03); v1.push_back(.05);
    v1.push_back(.01);
    v2.push_back(7);
    v4.push_back(0);

    v1.push_back(.125); v1.push_back(.175); v1.push_back(.3);
    v1.push_back(.19); v1.push_back(.21);
    v2.push_back(5);
    v4.push_back(0);

    v1.push_back(.1); v1.push_back(.2); v1.push_back(.01); v1.push_back(.02);
    v1.push_back(.03); v1.push_back(.04); v1.push_back(.6);
    v2.push_back(7);
    v4.push_back(0);

    v3.push_back(4);

    // VJ ins len
    v1.push_back(.05); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2); v1.push_back(.25); v1.push_back(.24); v1.push_back(.1);
    v2.push_back(7);
    v4.push_back(0);

    v3.push_back(7);

    // VJ ins nuc
    // prev A
    v1.push_back(.05); v1.push_back(.08); v1.push_back(.03); v1.push_back(.84);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(8);

    // prev C
    v1.push_back(.4); v1.push_back(.1); v1.push_back(.3); v1.push_back(.2);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(9);

    // prev G
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(10);

    // prev T
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.25); v1.push_back(.3);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(11);

    ModelParameterVector mvec(VJ_RECOMB, v1, v2, v3, v4);

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

    ClonotypeBuilder cl_builder;
    // CCCG.AC.GGTTT
    // poses:
    // vs: 0-1-2-3-4-5
    // js: 7-8-9-10-11
    cl_builder.setSequence("CCCGACGGTTT")
            .setNucleotideSeq()
            .addValignment(1, 4)
            .addValignment(3, 5)
            .addJalignment(1, 8)
            .addJalignment(2, 9)
            .addJalignment(3, 7);
    Clonotype clonotype = cl_builder.buildClonotype();

    MAAG maag = maag_builder.build(clonotype, true);

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

    YMIR_ASSERT2(maag.event_index(3, 2, 0, 0), mvec.event_index(VJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_index(3, 2, 1, 0), mvec.event_index(VJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(3, 2, 2, 0), mvec.event_index(VJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(3, 2, 3, 0), mvec.event_index(VJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(3, 2, 4, 0), mvec.event_index(VJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(3, 2, 5, 0), mvec.event_index(VJ_JOI_DEL, 2, 6))

    YMIR_ASSERT2(maag.event_probability(0, 0, 0, 0), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 0, 0))
    YMIR_ASSERT2(maag.event_probability(0, 0, 0, 1), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 0, 1))
    YMIR_ASSERT2(maag.event_probability(0, 0, 1, 0), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 2, 0))
    YMIR_ASSERT2(maag.event_probability(0, 0, 1, 2), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 2, 2))

    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 0), mvec.event_prob(VJ_VAR_DEL, 0, 4));
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 1), mvec.event_prob(VJ_VAR_DEL, 0, 3));
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 2), mvec.event_prob(VJ_VAR_DEL, 0, 2));
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 3), mvec.event_prob(VJ_VAR_DEL, 0, 1));
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 4), mvec.event_prob(VJ_VAR_DEL, 0, 0));
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 5), 0)

    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 0), mvec.event_prob(VJ_VAR_DEL, 2, 6))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 1), mvec.event_prob(VJ_VAR_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 2), mvec.event_prob(VJ_VAR_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 3), mvec.event_prob(VJ_VAR_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 4), mvec.event_prob(VJ_VAR_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 5), mvec.event_prob(VJ_VAR_DEL, 2, 1))

    YMIR_ASSERT2(maag.event_probability(2, 0, 0, 1), 0)

    YMIR_ASSERT2(maag.event_probability(3, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.event_probability(3, 0, 1, 0), mvec.event_prob(VJ_JOI_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_probability(3, 0, 2, 0), mvec.event_prob(VJ_JOI_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_probability(3, 0, 3, 0), mvec.event_prob(VJ_JOI_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_probability(3, 0, 4, 0), mvec.event_prob(VJ_JOI_DEL, 0, 5))
    YMIR_ASSERT2(maag.event_probability(3, 0, 5, 0), mvec.event_prob(VJ_JOI_DEL, 0, 6))

    YMIR_ASSERT2(maag.event_probability(3, 2, 0, 0), mvec.event_prob(VJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_probability(3, 2, 1, 0), mvec.event_prob(VJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_probability(3, 2, 2, 0), mvec.event_prob(VJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_probability(3, 2, 3, 0), mvec.event_prob(VJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_probability(3, 2, 4, 0), mvec.event_prob(VJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_probability(3, 2, 5, 0), mvec.event_prob(VJ_JOI_DEL, 2, 6))
YMIR_TEST_END


YMIR_TEST_START(test_maag_vdj)
    vector<prob_t> v1;  // param vec
    vector<eventind_t> v2;  // lens vec
    vector<eventind_t> v3;  // event classes
    vector<seq_len_t> v4;  // event family col numbers

    // V
    v1.push_back(.5); v1.push_back(.25); v1.push_back(.25);
    v2.push_back(3);

    v3.push_back(0);
    v4.push_back(0);

    // J-D
    // J-D (3 Js - 3 Ds)
    v1.push_back(.01); v1.push_back(.02); v1.push_back(.03);
    v1.push_back(.04); v1.push_back(.05); v1.push_back(.07);
    v1.push_back(.08); v1.push_back(.09); v1.push_back(.61);

    v2.push_back(9);

    v3.push_back(1);
    v4.push_back(3);

    // V del
    v1.push_back(.4); v1.push_back(.5); v1.push_back(.05); v1.push_back(.02); v1.push_back(.03);
    v2.push_back(5);

    v1.push_back(.3); v1.push_back(.1); v1.push_back(.2); v1.push_back(.4);
    v2.push_back(4);

    v1.push_back(.75); v1.push_back(.005); v1.push_back(.01); v1.push_back(.02);
    v1.push_back(.03); v1.push_back(.04); v1.push_back(.145);
    v2.push_back(7);

    v3.push_back(2);
    v4.push_back(0);
    v4.push_back(0);
    v4.push_back(0);

    // J del
    v1.push_back(.07); v1.push_back(.2); v1.push_back(.3);
    v1.push_back(.34); v1.push_back(.03); v1.push_back(.05);
    v1.push_back(.01);
    v2.push_back(7);

    v1.push_back(.125); v1.push_back(.175); v1.push_back(.3);
    v1.push_back(.19); v1.push_back(.21);
    v2.push_back(5);

    v1.push_back(.1); v1.push_back(.2); v1.push_back(.01); v1.push_back(.02);
    v1.push_back(.03); v1.push_back(.04); v1.push_back(.6);
    v2.push_back(7);

    v3.push_back(5);
    v4.push_back(0);
    v4.push_back(0);
    v4.push_back(0);

    // D1 dels
    v1.push_back(.17); // first row
    v1.push_back(.27);
    v1.push_back(.37); // second row
    v1.push_back(.19);

    v2.push_back(4);
    v4.push_back(2);

    // D2 dels
    // 3 rows 2 columns
    v1.push_back(.11); v1.push_back(.12);
    v1.push_back(.13); v1.push_back(.14);
    v1.push_back(.15); v1.push_back(.35);

    v2.push_back(6);
    v4.push_back(2);

    // D3 dels
    v1.push_back(.11); // first row
    v1.push_back(.22);
    v1.push_back(.33); // second row
    v1.push_back(.34);

    v2.push_back(4);
    v4.push_back(2);

    v3.push_back(8);

    // VD ins len
    v1.push_back(.05); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2); v1.push_back(.25); v1.push_back(.24); v1.push_back(.1);
    v2.push_back(7);

    v3.push_back(11);
    v4.push_back(0);

    // DJ ins len
    v1.push_back(.1); v1.push_back(.24); v1.push_back(.25); v1.push_back(.05); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2);
    v2.push_back(7);

    v3.push_back(12);
    v4.push_back(0);

    // VD ins nuc
    // prev A
    v1.push_back(.05); v1.push_back(.08); v1.push_back(.03); v1.push_back(.84);
    v2.push_back(4);

    v3.push_back(13);
    v4.push_back(0);

    // prev C
    v1.push_back(.4); v1.push_back(.1); v1.push_back(.3); v1.push_back(.2);
    v2.push_back(4);

    v3.push_back(14);
    v4.push_back(0);

    // prev G
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2);
    v2.push_back(4);

    v3.push_back(15);
    v4.push_back(0);

    // prev T
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.25); v1.push_back(.3);
    v2.push_back(4);

    v3.push_back(16);
    v4.push_back(0);

    // DJ ins nuc
    // prev A
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2);
    v2.push_back(4);

    v3.push_back(17);
    v4.push_back(0);

    // prev C
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.25); v1.push_back(.3);
    v2.push_back(4);

    v3.push_back(18);
    v4.push_back(0);

    // prev G
    v1.push_back(.05); v1.push_back(.08); v1.push_back(.03); v1.push_back(.84);
    v2.push_back(4);

    v3.push_back(19);
    v4.push_back(0);

    // prev T
    v1.push_back(.4); v1.push_back(.1); v1.push_back(.3); v1.push_back(.2);
    v2.push_back(4);

    v3.push_back(20);
    v4.push_back(0);

    ModelParameterVector mvec(VDJ_RECOMB, v1, v2, v3, v3);

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

    ClonotypeBuilder cl_builder;
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
            .setNucleotideSeq()
            .addValignment(1, 4)
            .addValignment(3, 5)
            .addJalignment(1, 8)
            .addJalignment(2, 9)
            .addJalignment(3, 7)
            .addDalignment(2, 2, 4, 2, 4)
            .addDalignment(2, 3, 6, 6, 9)
            .addDalignment(3, 5, 7, 4, 6)
            .addDalignment(3, 1, 4, 1, 4)
            .addDalignment(3, 3, 5, 6, 8)
            .addDalignment(1, 1, 4, 8, 11);
    Clonotype clonotype = cl_builder.buildClonotype();

    MAAG maag = maag_builder.build(clonotype, true);

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

    YMIR_ASSERT2(maag.event_index(2, 0, 0, 10), 0)
    YMIR_ASSERT2(maag.event_index(2, 0, 5, 10), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 5))
    YMIR_ASSERT2(maag.event_index(2, 0, 0, 0), mvec.event_index(VDJ_VAR_DIV_INS_LEN, 0, 0))

    YMIR_ASSERT2(maag.event_index(3, 0, 1, 3), mvec.event_index(VDJ_DIV_DEL, 1, 1, 2))
    YMIR_ASSERT2(maag.event_index(3, 0, 5, 7), mvec.event_index(VDJ_DIV_DEL, 1, 2, 1))
    YMIR_ASSERT2(maag.event_index(3, 0, 6, 8), mvec.event_index(VDJ_DIV_DEL, 1, 3, 0))
    YMIR_ASSERT2(maag.event_index(3, 0, 5, 8), mvec.event_index(VDJ_DIV_DEL, 1, 2, 0))

    YMIR_ASSERT2(maag.event_index(3, 2, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 2, 5, 0), 0)
    YMIR_ASSERT2(maag.event_index(3, 2, 0, 3), 0)
    YMIR_ASSERT2(maag.event_index(3, 2, 7, 10), mvec.event_index(VDJ_DIV_DEL, 0, 0, 0))
    YMIR_ASSERT2(maag.event_index(3, 2, 8, 10), mvec.event_index(VDJ_DIV_DEL, 0, 1, 0))
    YMIR_ASSERT2(maag.event_index(3, 2, 7, 9), mvec.event_index(VDJ_DIV_DEL, 0, 0, 1))

    YMIR_ASSERT2(maag.event_index(4, 0, 0, 0), mvec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 5))
    YMIR_ASSERT2(maag.event_index(4, 0, 5, 4), mvec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 4))
    YMIR_ASSERT2(maag.event_index(4, 0, 5, 5), 0)
    YMIR_ASSERT2(maag.event_index(4, 0, 6, 1), mvec.event_index(VDJ_DIV_JOI_INS_LEN, 0, 0))
    YMIR_ASSERT2(maag.event_index(4, 0, 10, 2), 0)

    YMIR_ASSERT2(maag.event_index(5, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.event_index(5, 0, 1, 0), mvec.event_index(VDJ_JOI_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_index(5, 0, 2, 0), mvec.event_index(VDJ_JOI_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_index(5, 0, 3, 0), mvec.event_index(VDJ_JOI_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_index(5, 0, 4, 0), mvec.event_index(VDJ_JOI_DEL, 0, 5))
    YMIR_ASSERT2(maag.event_index(5, 0, 5, 0), mvec.event_index(VDJ_JOI_DEL, 0, 6))

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


    YMIR_ASSERT2(maag.event_probability(0, 0, 0, 0), mvec.event_prob(VDJ_VAR_GEN, 0, 0))
    YMIR_ASSERT2(maag.event_probability(0, 1, 0, 0), mvec.event_prob(VDJ_VAR_GEN, 0, 2))

    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 0), mvec.event_prob(VDJ_VAR_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 1), mvec.event_prob(VDJ_VAR_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 2), mvec.event_prob(VDJ_VAR_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 3), mvec.event_prob(VDJ_VAR_DEL, 0, 1))
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 4), mvec.event_prob(VDJ_VAR_DEL, 0, 0))
    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 5), 0)

    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 0), mvec.event_prob(VDJ_VAR_DEL, 2, 6))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 1), mvec.event_prob(VDJ_VAR_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 2), mvec.event_prob(VDJ_VAR_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 3), mvec.event_prob(VDJ_VAR_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 4), mvec.event_prob(VDJ_VAR_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 5), mvec.event_prob(VDJ_VAR_DEL, 2, 1))

    YMIR_ASSERT2(maag.event_probability(2, 0, 0, 10), 0)

    YMIR_ASSERT2(maag.event_probability(3, 0, 1, 3), mvec.event_prob(VDJ_DIV_DEL, 1, 1, 2))
    YMIR_ASSERT2(maag.event_probability(3, 0, 5, 7), mvec.event_prob(VDJ_DIV_DEL, 1, 2, 1))
    YMIR_ASSERT2(maag.event_probability(3, 0, 6, 8), mvec.event_prob(VDJ_DIV_DEL, 1, 3, 0))
    YMIR_ASSERT2(maag.event_probability(3, 0, 5, 8), mvec.event_prob(VDJ_DIV_DEL, 1, 2, 0))

    YMIR_ASSERT2(maag.event_probability(3, 2, 0, 0), 0)
    YMIR_ASSERT2(maag.event_probability(3, 2, 5, 0), 0)
    YMIR_ASSERT2(maag.event_probability(3, 2, 0, 3), 0)

    YMIR_ASSERT2(maag.event_probability(4, 0, 5, 5), 0)
    YMIR_ASSERT2(maag.event_probability(4, 0, 10, 2), 0)

    YMIR_ASSERT2(maag.event_probability(5, 0, 0, 0), 0)
    YMIR_ASSERT2(maag.event_probability(5, 0, 1, 0), mvec.event_prob(VDJ_JOI_DEL, 0, 2))
    YMIR_ASSERT2(maag.event_probability(5, 0, 2, 0), mvec.event_prob(VDJ_JOI_DEL, 0, 3))
    YMIR_ASSERT2(maag.event_probability(5, 0, 3, 0), mvec.event_prob(VDJ_JOI_DEL, 0, 4))
    YMIR_ASSERT2(maag.event_probability(5, 0, 4, 0), mvec.event_prob(VDJ_JOI_DEL, 0, 5))
    YMIR_ASSERT2(maag.event_probability(5, 0, 5, 0), mvec.event_prob(VDJ_JOI_DEL, 0, 6))

    YMIR_ASSERT2(maag.event_probability(5, 2, 0, 0), mvec.event_prob(VDJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_probability(5, 2, 1, 0), mvec.event_prob(VDJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_probability(5, 2, 2, 0), mvec.event_prob(VDJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_probability(5, 2, 3, 0), mvec.event_prob(VDJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_probability(5, 2, 4, 0), mvec.event_prob(VDJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_probability(5, 2, 5, 0), mvec.event_prob(VDJ_JOI_DEL, 2, 6))

    YMIR_ASSERT2(maag.event_probability(6, 0, 0, 0), mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 0, 1))
    YMIR_ASSERT2(maag.event_probability(6, 0, 0, 1), mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 0, 2))
    YMIR_ASSERT2(maag.event_probability(6, 0, 0, 2), mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 0, 0))
    YMIR_ASSERT2(maag.event_probability(6, 0, 2, 0), mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 2, 1))
    YMIR_ASSERT2(maag.event_probability(6, 0, 2, 1), mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 2, 2))
    YMIR_ASSERT2(maag.event_probability(6, 0, 2, 2), mvec.event_prob(VDJ_JOI_DIV_GEN, 0, 2, 0))

    // i don't want to compute by hand all this crazy matrices!
    // i've already tested chain products in previous tests!
    // ):<
    // also it's (A, B) not (A - B < eps) because this results are pretty precise on this toy example
    YMIR_ASSERT2(maag.fullProbability(0, 0, 0), maag_builder.build(clonotype, false).fullProbability(0, 0, 0))
    YMIR_ASSERT2(maag.fullProbability(1, 1, 1), maag_builder.build(clonotype, false).fullProbability(1, 1, 1))
    YMIR_ASSERT2(maag.fullProbability(0, 2, 2), maag_builder.build(clonotype, false).fullProbability(0, 2, 2))

YMIR_TEST_END


YMIR_TEST_START(test_maag_builder_replace)
    YMIR_ASSERT(false)
YMIR_TEST_END


YMIR_TEST_START(test_model_vj_file)
    ProbabilisticAssemblingModel model1(TEST_DATA_FOLDER + "randomfile");
    YMIR_ASSERT(!model1.status())

    vector<prob_t> v1;  // param vec
    vector<eventind_t> v2;  // lens vec
    vector<eventind_t> v3;  // event classes
    vector<seq_len_t> v4;  // event family col numbers

    // V-J
    v1.push_back(.05); v1.push_back(.025); v1.push_back(.035); // J1
    v1.push_back(.045); v1.push_back(.055); v1.push_back(.065); // J2
    v1.push_back(.075); v1.push_back(.085); v1.push_back(.565); // J3
    v2.push_back(9);

    v3.push_back(0);
    v4.push_back(3);

    // V del
    v1.push_back(.4); v1.push_back(.5); v1.push_back(.05); v1.push_back(.02); v1.push_back(.03);
    v2.push_back(5);

    v1.push_back(.3); v1.push_back(.1); v1.push_back(.2); v1.push_back(.4);
    v2.push_back(4);

    v1.push_back(.75); v1.push_back(.005); v1.push_back(.01); v1.push_back(.02);
    v1.push_back(.03); v1.push_back(.04); v1.push_back(.145);
    v2.push_back(7);

    v3.push_back(1);
    v4.push_back(0);
    v4.push_back(0);
    v4.push_back(0);

    // J del
    v1.push_back(.07); v1.push_back(.2); v1.push_back(.3);
    v1.push_back(.34); v1.push_back(.03); v1.push_back(.05);
    v1.push_back(.01);
    v2.push_back(7);
    v4.push_back(0);

    v1.push_back(.125); v1.push_back(.175); v1.push_back(.3);
    v1.push_back(.19); v1.push_back(.21);
    v2.push_back(5);
    v4.push_back(0);

    v1.push_back(.1); v1.push_back(.2); v1.push_back(.01); v1.push_back(.02);
    v1.push_back(.03); v1.push_back(.04); v1.push_back(.6);
    v2.push_back(7);
    v4.push_back(0);

    v3.push_back(4);

    // VJ ins len
    v1.push_back(.05); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2); v1.push_back(.25); v1.push_back(.24); v1.push_back(.01);
    v2.push_back(7);
    v4.push_back(0);

    v3.push_back(7);

    // VJ ins nuc
    // prev A
    v1.push_back(.05); v1.push_back(.08); v1.push_back(.03); v1.push_back(.84);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(8);

    // prev C
    v1.push_back(.4); v1.push_back(.1); v1.push_back(.3); v1.push_back(.2);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(9);

    // prev G
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.15); v1.push_back(.5);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(10);

    // prev T
    v1.push_back(.15); v1.push_back(.1); v1.push_back(.25); v1.push_back(.5);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(11);

    ModelParameterVector mvec(VJ_RECOMB, v1, v2, v3, v4);

    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vj_model/");
    YMIR_ASSERT(model.status())

    YMIR_ASSERT(mvec == model.event_probabilities())


YMIR_TEST_END


YMIR_TEST_START(test_model_vdj_file)
    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vdj_model/");
    YMIR_ASSERT(model.status())

YMIR_TEST_END


YMIR_TEST_START(test_model_vj_maag)
    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vj_model/");
    YMIR_ASSERT(model.status())

//    MAAG maag = model.buildGraphs(cloneset, true, false)[1];
//
//    YMIR_ASSERT2(maag.event_index(0, 0, 0, 0), mvec.event_index(VJ_VAR_JOI_GEN, 0, 0, 0))
//    YMIR_ASSERT2(maag.event_index(0, 0, 0, 1), mvec.event_index(VJ_VAR_JOI_GEN, 0, 0, 1))
//    YMIR_ASSERT2(maag.event_index(0, 0, 1, 0), mvec.event_index(VJ_VAR_JOI_GEN, 0, 2, 0))
//    YMIR_ASSERT2(maag.event_index(0, 0, 1, 2), mvec.event_index(VJ_VAR_JOI_GEN, 0, 2, 2))
//
//    YMIR_ASSERT2(maag.event_index(1, 0, 0, 0), mvec.event_index(VJ_VAR_DEL, 0, 4));
//    YMIR_ASSERT2(maag.event_index(1, 0, 0, 1), mvec.event_index(VJ_VAR_DEL, 0, 3));
//    YMIR_ASSERT2(maag.event_index(1, 0, 0, 2), mvec.event_index(VJ_VAR_DEL, 0, 2));
//    YMIR_ASSERT2(maag.event_index(1, 0, 0, 3), mvec.event_index(VJ_VAR_DEL, 0, 1));
//    YMIR_ASSERT2(maag.event_index(1, 0, 0, 4), mvec.event_index(VJ_VAR_DEL, 0, 0));
//    YMIR_ASSERT2(maag.event_index(1, 0, 0, 5), 0)
//
//    YMIR_ASSERT2(maag.event_index(1, 1, 0, 0), mvec.event_index(VJ_VAR_DEL, 2, 6))
//    YMIR_ASSERT2(maag.event_index(1, 1, 0, 1), mvec.event_index(VJ_VAR_DEL, 2, 5))
//    YMIR_ASSERT2(maag.event_index(1, 1, 0, 2), mvec.event_index(VJ_VAR_DEL, 2, 4))
//    YMIR_ASSERT2(maag.event_index(1, 1, 0, 3), mvec.event_index(VJ_VAR_DEL, 2, 3))
//    YMIR_ASSERT2(maag.event_index(1, 1, 0, 4), mvec.event_index(VJ_VAR_DEL, 2, 2))
//    YMIR_ASSERT2(maag.event_index(1, 1, 0, 5), mvec.event_index(VJ_VAR_DEL, 2, 1))
//
//    YMIR_ASSERT2(maag.event_index(2, 0, 0, 0), mvec.event_index(VJ_VAR_JOI_INS_LEN, 0, maag.position(6) - maag.position(0) - 1))
//    YMIR_ASSERT2(maag.event_index(2, 0, 0, 1), 0)
//
//    YMIR_ASSERT2(maag.event_index(3, 0, 0, 0), 0)
//    YMIR_ASSERT2(maag.event_index(3, 0, 1, 0), mvec.event_index(VJ_JOI_DEL, 0, 2))
//    YMIR_ASSERT2(maag.event_index(3, 0, 2, 0), mvec.event_index(VJ_JOI_DEL, 0, 3))
//    YMIR_ASSERT2(maag.event_index(3, 0, 3, 0), mvec.event_index(VJ_JOI_DEL, 0, 4))
//    YMIR_ASSERT2(maag.event_index(3, 0, 4, 0), mvec.event_index(VJ_JOI_DEL, 0, 5))
//    YMIR_ASSERT2(maag.event_index(3, 0, 5, 0), mvec.event_index(VJ_JOI_DEL, 0, 6))
//
//    YMIR_ASSERT2(maag.event_index(3, 2, 0, 0), mvec.event_index(VJ_JOI_DEL, 2, 1))
//    YMIR_ASSERT2(maag.event_index(3, 2, 1, 0), mvec.event_index(VJ_JOI_DEL, 2, 2))
//    YMIR_ASSERT2(maag.event_index(3, 2, 2, 0), mvec.event_index(VJ_JOI_DEL, 2, 3))
//    YMIR_ASSERT2(maag.event_index(3, 2, 3, 0), mvec.event_index(VJ_JOI_DEL, 2, 4))
//    YMIR_ASSERT2(maag.event_index(3, 2, 4, 0), mvec.event_index(VJ_JOI_DEL, 2, 5))
//    YMIR_ASSERT2(maag.event_index(3, 2, 5, 0), mvec.event_index(VJ_JOI_DEL, 2, 6))
//
//    YMIR_ASSERT2(maag.event_probability(0, 0, 0, 0), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 0, 0))
//    YMIR_ASSERT2(maag.event_probability(0, 0, 0, 1), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 0, 1))
//    YMIR_ASSERT2(maag.event_probability(0, 0, 1, 0), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 2, 0))
//    YMIR_ASSERT2(maag.event_probability(0, 0, 1, 2), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 2, 2))
//
//    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 0), mvec.event_prob(VJ_VAR_DEL, 0, 4));
//    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 1), mvec.event_prob(VJ_VAR_DEL, 0, 3));
//    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 2), mvec.event_prob(VJ_VAR_DEL, 0, 2));
//    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 3), mvec.event_prob(VJ_VAR_DEL, 0, 1));
//    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 4), mvec.event_prob(VJ_VAR_DEL, 0, 0));
//    YMIR_ASSERT2(maag.event_probability(1, 0, 0, 5), 0)
//
//    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 0), mvec.event_prob(VJ_VAR_DEL, 2, 6))
//    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 1), mvec.event_prob(VJ_VAR_DEL, 2, 5))
//    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 2), mvec.event_prob(VJ_VAR_DEL, 2, 4))
//    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 3), mvec.event_prob(VJ_VAR_DEL, 2, 3))
//    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 4), mvec.event_prob(VJ_VAR_DEL, 2, 2))
//    YMIR_ASSERT2(maag.event_probability(1, 1, 0, 5), mvec.event_prob(VJ_VAR_DEL, 2, 1))
//
//    YMIR_ASSERT2(maag.event_probability(2, 0, 0, 1), 0)
//
//    YMIR_ASSERT2(maag.event_probability(3, 0, 0, 0), 0)
//    YMIR_ASSERT2(maag.event_probability(3, 0, 1, 0), mvec.event_prob(VJ_JOI_DEL, 0, 2))
//    YMIR_ASSERT2(maag.event_probability(3, 0, 2, 0), mvec.event_prob(VJ_JOI_DEL, 0, 3))
//    YMIR_ASSERT2(maag.event_probability(3, 0, 3, 0), mvec.event_prob(VJ_JOI_DEL, 0, 4))
//    YMIR_ASSERT2(maag.event_probability(3, 0, 4, 0), mvec.event_prob(VJ_JOI_DEL, 0, 5))
//    YMIR_ASSERT2(maag.event_probability(3, 0, 5, 0), mvec.event_prob(VJ_JOI_DEL, 0, 6))
//
//    YMIR_ASSERT2(maag.event_probability(3, 2, 0, 0), mvec.event_prob(VJ_JOI_DEL, 2, 1))
//    YMIR_ASSERT2(maag.event_probability(3, 2, 1, 0), mvec.event_prob(VJ_JOI_DEL, 2, 2))
//    YMIR_ASSERT2(maag.event_probability(3, 2, 2, 0), mvec.event_prob(VJ_JOI_DEL, 2, 3))
//    YMIR_ASSERT2(maag.event_probability(3, 2, 3, 0), mvec.event_prob(VJ_JOI_DEL, 2, 4))
//    YMIR_ASSERT2(maag.event_probability(3, 2, 4, 0), mvec.event_prob(VJ_JOI_DEL, 2, 5))
//    YMIR_ASSERT2(maag.event_probability(3, 2, 5, 0), mvec.event_prob(VJ_JOI_DEL, 2, 6))

YMIR_TEST_END


YMIR_TEST_START(test_model_vdj_maag)
    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vdj_model/");
    YMIR_ASSERT(model.status())

YMIR_TEST_END


struct TestInfo {
    string test_name;
    vector<string> failed_cases;

    TestInfo(const string& name, const vector<string> vec) : test_name(name), failed_cases(vec) {}
};


int main() {
    //
    // MEGA TO-DO: make good tests with some unit-testing framework (CTest)
    //

    //**************  INITIALISATION  **************//
    size_t tests_passed = 0, all_tests = 0;
    vector<TestInfo> failed_test_info;
    //**************  **************//



    //**************  TEST CASES  **************//
    YMIR_TEST(test_basic(), "basic test")

    // Tests for ModelParameterVector
    YMIR_TEST(test_model_param_vec_vj(), "ModelParameterVector VJ access w/o laplace")
    YMIR_TEST(test_model_param_vec_vdj(), "ModelParameterVector VDJ access w/o laplace")

    // Tests for gene segments classes
    YMIR_TEST(test_genesegmentalphabet(), "GeneSegmentAlphabet initialisation and access")
    YMIR_TEST(test_vdjgenes2(), "VDJRecombinationGenes, 2 gene segments, initialisation and access")
    YMIR_TEST(test_vdjgenes3(), "VDJRecombinationGenes, 3 gene segments, initialisation and access")
    YMIR_TEST(test_genesegmentalphabet_read(), "GeneSegmentAlphabet file reading / writing")
    YMIR_TEST(test_vdjgenes_read(), "VDJRecombinationGenes file reading / writing")

    // Tests for clone, clone alignment and clone builder classes.
    YMIR_TEST(test_clone(), "Clonotype, constructor / access")
    YMIR_TEST(test_clonebuilder_clonealign(), "ClonotypeBuilder, Clonotype building")

    // Tests for default naive sequences aligners.
    YMIR_TEST(test_nuc_aligner(), "Naive nucleotide sequence aligner")
    YMIR_TEST(test_aa_aligner(), "Naive amino acid sequence aligner")

    // Test for MiTCR parser.
    YMIR_TEST(test_mitcr_vj(), "MiTCR parser test for alpha chain")
    YMIR_TEST(test_mitcr_vdj(), "MiTCR parser test for beta chain")

    // Tests for clonal repertoires and clonal repertoire views.
    YMIR_TEST(test_clorep(), "Cloneset creating / access")
    YMIR_TEST(test_clorep_view(), "ClonesetView creating / access")

    // Tests for markov chain.
    YMIR_TEST(test_markovchain_nuc(), "Markov chain nucleotide")
    YMIR_TEST(test_markovchain_aa(), "Markov chain amino acid")

    // Test for Multi-Matrix Chains
    YMIR_TEST(test_mmc(), "Multi-Matrix chain all interface")

//    // Tests for MAAG / MAAG builder
    YMIR_TEST(test_maag_vj(), "MAAG VJ building and computing")
    YMIR_TEST(test_maag_vdj(), "MAAG VDJ building and computing")
    YMIR_TEST(test_maag_builder_replace(), "MAAG Builder replace event probabilities")

    // Tests for assembling statistical model (ASM) reading / writing files.
    YMIR_TEST(test_model_vj_file(), "VJ Model constructing from a file")
//    YMIR_TEST(test_model_vdj_file(), "VDJ Model constructing from a file")
//    YMIR_TEST(test_model_vj_maag(), "VJ Model creating MAAGs")
//    YMIR_TEST(test_model_vdj_maag(), "VDJ Model creating MAAGs")

    // Test for computing full nucleotide probabilities of repertoire with ASM.

    // Test for computing full amino acid probabilities of repertoire with ASM.

    // Tests for statistical inference of ASM parameters.

    //**************  **************//



    //**************  TESTING RESULTS  **************//
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