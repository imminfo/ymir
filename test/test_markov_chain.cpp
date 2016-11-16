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


YMIR_TEST_START(test_markovchain_nuc_mono)

    vector<prob_t> probs = {.1, .2, .3, .4};
    MonoNucInsertionModel m(probs.begin());

    string s = "ACGT";

    // .1 * .2 * .3 * .4
    YMIR_ASSERT(abs(m.nucProbability(s.begin(), 4, NULL_CHAR)) - .0024 < 1e-18)
    YMIR_ASSERT(abs(m.nucProbability(s.begin(), 4)) - .0024 < 1e-18)
    YMIR_ASSERT(abs(m.nucProbability(s, NULL_CHAR)) - .0024 < 1e-18)

    // .1 * .2 * .3 * .4
    YMIR_ASSERT(abs(m.nucProbability(s.begin(), 4, 'A')) - .0024 < 1e-18)
    YMIR_ASSERT(abs(m.nucProbability(s, 'A')) - .0024 < 1e-18)

    // 1 * .2 * .3
    YMIR_ASSERT(abs(m.nucProbability(s.substr(0, 3), 'A')) - .006 < 1e-18)
    // .1 * .2 * .3
    YMIR_ASSERT(abs(m.nucProbability(s.substr(0, 3), NULL_CHAR)) - .006 < 1e-18)

    probs = {1, 0, 0, 0};
    m = MonoNucInsertionModel(probs.begin());
    std::default_random_engine rg;
    YMIR_ASSERT2(m.generate(5, rg), "AAAAA");

    probs = {0, 0, 0, 1};
    m = MonoNucInsertionModel(probs.begin());
    YMIR_ASSERT2(m.generate(1, rg), "T");

    probs = {0, 1, 0, 0};
    m = MonoNucInsertionModel(probs.begin());
    YMIR_ASSERT2(m.generate(3, rg), "CCC");

YMIR_TEST_END


YMIR_TEST_START(test_markovchain_nuc_di)

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

    DiNucInsertionModel mc(mat);
    string s = "ACGT";

    // .25 * .7 * .3 * .2 = .0105
    YMIR_ASSERT(mc.nucProbability("") == 1);
    YMIR_ASSERT(mc.nucProbability(s) - .0105 < 1e-16);
    YMIR_ASSERT(mc.nucProbability(s.begin(), 4) - .0105 < 1e-16);

    // .4 * .7 * .3 * .2 = .0126
    YMIR_ASSERT(mc.nucProbability(s.begin(), 4, 'C') - .0168 < 1e-16);

    mc = DiNucInsertionModel(vec.begin());

    YMIR_ASSERT(mc.nucProbability(s) - .0105 < 1e-18);
    YMIR_ASSERT(mc.nucProbability(s.begin(), 4, 'C') - .0168 < 1e-16);

    // .25 * .25 * .1
    YMIR_ASSERT(mc.nucProbability(s.rbegin(), 3, '_') - .00625 < 1e-16);

YMIR_TEST_END


YMIR_TEST_START(test_markovchain_aa_mono)

    vector<prob_t> probs = {.1, .2, .3, .4};
    MonoNucInsertionModel m(probs.begin(), .5);

    //
    // exact same aminoacid
    //
    // {'M', "ATG"}
    YMIR_ASSERT3(m.aaProbability("M", 1, 1, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("M", 1, 2, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("M", 1, 3, 0, 0), 0)

    YMIR_ASSERT3(m.aaProbability("M", 1, 1, 63, 63), .1)
    YMIR_ASSERT3(m.aaProbability("M", 1, 2, 63, 63), .1 * .4)
    YMIR_ASSERT3(m.aaProbability("M", 1, 3, 63, 63), .1 * .4 * .3)
    YMIR_ASSERT3(m.aaProbability("M", 2, 3, 63, 63), .4 * .3)
    YMIR_ASSERT3(m.aaProbability("M", 3, 3, 63, 63), .3)

    //
    // {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 0, 0), 0)

    // 111000
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 56, 56), 3 * .1)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 56, 56), .1*.4 + .1*.4 + .1*.4)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 56, 56), .1*.4*.4 + .1*.4*.2 + .1*.4*.1)
    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 56, 56), (3 * .4))
    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 56, 56), .4 * (.4 + .2 + .1))
    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 56, 56), (.1 + .2 + .4))

    // 101000
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 40, 40), 2 * .1)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 40, 40), .1*.4 + .1*.4)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 40, 40), .1*.4*.4 + .1*.4*.1)
    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 40, 40), .4*2)
    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 40, 40), .4*.4 + .4*.1)
    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 40, 40), .4 + .1)

    // 110000 and 011000 -> 010000
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 48, 24), .1)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 48, 24), .1 * .4)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 48, 24), .1 * .4 * .2)
    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 48, 24), .4)
    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 48, 24), .4 * .2)
    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 48, 24), .2)


    //
    // neighbour aminoacids
    //
    // {'A', "GCT"}, {'A', "GCC"}, {'A', "GCA"}, {'A', "GCG"}
    // {'F', "TTT"}, {'F', "TTC"}
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("AF", 1, 5, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("AF", 1, 6, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("AF", 2, 4, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("AF", 2, 5, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("AF", 2, 6, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("AF", 3, 4, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("AF", 3, 5, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("AF", 3, 6, 0, 0), 0)

    // 111100 and 110000 (full)
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 60, 48), 0.048)
    YMIR_ASSERT3(m.aaProbability("AF", 1, 5, 60, 48), 0.0192)
    YMIR_ASSERT3(m.aaProbability("AF", 1, 6, 60, 48), 0.00576)
    YMIR_ASSERT3(m.aaProbability("AF", 2, 4, 60, 48), 0.16)
    YMIR_ASSERT3(m.aaProbability("AF", 2, 5, 60, 48), 0.064)
    YMIR_ASSERT3(m.aaProbability("AF", 2, 6, 60, 48), 0.0192)
    YMIR_ASSERT3(m.aaProbability("AF", 3, 4, 60, 48), .8)
    YMIR_ASSERT3(m.aaProbability("AF", 3, 5, 60, 48), 0.32)
    YMIR_ASSERT3(m.aaProbability("AF", 3, 6, 60, 48), 0.096)


    //
    // distant aminoacids
    //
    // {'N', "AAT"}, {'N', "AAC"}
    // {'G', "GGT"}, {'G', "GGC"}, {'G', "GGA"}, {'G', "GGG"}
    // {'Y', "TAT"}, {'Y', "TAC"}
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 7, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 8, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 9, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 7, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 8, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 9, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 7, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 8, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 9, 0, 0), 0)
    
    // 110000, 111100, 110000 - full
    // G - (4*.3) * (4*.3) * 1 = 1.44
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 7, 48, 48), 0.000432)
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 8, 48, 48), (.1*.1*.4 + .1*.1*.2) * (.3*.3) * (.4*.1 + .4*.1))
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 9, 48, 48), (.1*.1*.4 + .1*.1*.2) * (.3*.3) * (.4*.1*.4 + .4*.1*.2))
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 7, 48, 48), 0.00432)
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 8, 48, 48), 0.000432)
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 9, 48, 48), 0.0001296)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 7, 48, 48), 0.0432)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 8, 48, 48), 0.00432)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 9, 48, 48), 0.001296)

    // 010000, 111100, 110000 - full
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 7, 16, 48), 0.000144)
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 8, 16, 48), (.1*.1*.2) * (.3*.3) * (.4*.1 + .4*.1))
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 9, 16, 48), (.1*.1*.2) * (.3*.3) * (.4*.1*.4 + .4*.1*.2))
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 7, 16, 48), 0.00144)
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 8, 16, 48), 0.000144)
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 9, 16, 48), (.1*.2) * (.3*.3) * (.4*.1*.4 + .4*.1*.2))
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 7, 16, 48), 0.0144)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 8, 16, 48), 0.00144)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 9, 16, 48), 0.000432)

    // 110000, 111100, 100000 - full
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 7, 48, 32), 0.000216)
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 8, 48, 32), (.1*.1*.4 + .1*.1*.2) * (.3*.3) * (.4*.1))
    YMIR_ASSERT3(m.aaProbability("NGY", 1, 9, 48, 32), (.1*.1*.4 + .1*.1*.2) * (.3*.3) * (.4*.1*.4))
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 7, 48, 32), 0.00216)
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 8, 48, 32), 0.000216)
    YMIR_ASSERT3(m.aaProbability("NGY", 2, 9, 48, 32), (.1*.4 + .1*.2) * (.3*.3) * (.4*.1*.4))
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 7, 48, 32), 0.0216)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 8, 48, 32), 0.00216)
    YMIR_ASSERT3(m.aaProbability("NGY", 3, 9, 48, 32), 0.000864)

YMIR_TEST_END


YMIR_TEST_START(test_markovchain_aa_di)

    event_matrix_t mat;
    mat.resize(4, 4);
    // A
    mat(0, 0) = .1; // A->A
    mat(1, 0) = .2; // C->A
    mat(2, 0) = .3; // G->A
    mat(3, 0) = .4; // T->A
    // C
    mat(0, 1) = .5; // A->C
    mat(1, 1) = .1; // C->C
    mat(2, 1) = .3; // G->C
    mat(3, 1) = .1; // T->C
    // G
    mat(0, 2) = .1; // A->G
    mat(1, 2) = .5; // C->G
    mat(2, 2) = .1; // G->G
    mat(3, 2) = .3; // T->G
    // T
    mat(0, 3) = .4; // A->T
    mat(1, 3) = .3; // C->T
    mat(2, 3) = .2; // G->T
    mat(3, 3) = .1; // T->T

    // prev -> next
    std::map<std::string, prob_t> m_;
    m_["AA"] = mat(0, 0);
    m_["CA"] = mat(1, 0);
    m_["GA"] = mat(2, 0);
    m_["TA"] = mat(3, 0);

    m_["AC"] = mat(0, 1);
    m_["CC"] = mat(1, 1);
    m_["GC"] = mat(2, 1);
    m_["TC"] = mat(3, 1);

    m_["AG"] = mat(0, 2);
    m_["CG"] = mat(1, 2);
    m_["GG"] = mat(2, 2);
    m_["TG"] = mat(3, 2);

    m_["AT"] = mat(0, 3);
    m_["CT"] = mat(1, 3);
    m_["GT"] = mat(2, 3);
    m_["TT"] = mat(3, 3);

    auto comp_prob = [](const std::string &str, const  std::map<std::string, prob_t> &m_) {
        prob_t res = 1;
        for (int i = 0; i < str.size() - 1; ++i) {
            res *= m_.find(str.substr(i, 2))->second;
        }
        return res;
    };
    

    DiNucInsertionModel m(mat, .5);

    //
    // exact same aminoacid
    //
    // {'M', "ATG"}
//    YMIR_ASSERT3(m.aaProbability("M", 1, 1, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("M", 1, 2, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("M", 1, 3, 0, 0), 0)
//
//    YMIR_ASSERT3(m.aaProbability("M", 0, 1, 63, 63), .25)
//    YMIR_ASSERT3(m.aaProbability("M", 1, 1, 63, 63), .25)
//    YMIR_ASSERT3(m.aaProbability("M", 1, 2, 63, 63), .25 * .4)
//    YMIR_ASSERT3(m.aaProbability("M", 1, 3, 63, 63), .25 * .4 * .3)
//    YMIR_ASSERT3(m.aaProbability("M", 2, 3, 63, 63), .4 * .3)
//    YMIR_ASSERT3(m.aaProbability("M", 3, 3, 63, 63), .3)
//
//    // {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 0, 0), 0)
//
//    // 111000
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 56, 56), m_["GA"] + m_["GA"] + m_["GA"])
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 56, 56), comp_prob("GAT", m_) * 3)
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 56, 56), comp_prob("GATT", m_) + comp_prob("GATC", m_) + comp_prob("GATA", m_))
//    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 56, 56), m_["AT"] + m_["AT"] + m_["AT"])
//    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 56, 56), comp_prob("ATT", m_) + comp_prob("ATC", m_) + comp_prob("ATA", m_))
//    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 56, 56), m_["TT"] + m_["TC"] + m_["TA"])
//
//    // 101000
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 40, 40), m_["GA"] + m_["GA"])
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 40, 40), comp_prob("GAT", m_) * 2)
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 40, 40), comp_prob("GATT", m_) + comp_prob("GATA", m_))
//    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 40, 40), m_["AT"] + m_["AT"])
//    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 40, 40), comp_prob("ATT", m_) + comp_prob("ATA", m_))
//    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 40, 40), m_["TT"] + m_["TA"])
//
//    // 110000 and 011000 -> 010000
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 48, 24), m_["GA"])
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 48, 24), comp_prob("GAT", m_))
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 48, 24), comp_prob("GATC", m_))
//    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 48, 24), m_["AT"])
//    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 48, 24), comp_prob("ATC", m_))
//    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 48, 24), m_["TC"])
//
//    // 100000 and 101000 -> 100000
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 32, 40), m_["GA"])
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 48, 40), comp_prob("GAT", m_))
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 48, 40), comp_prob("GATT", m_))
//
//    // prev: {'H', "CAT"}, {'H', "CAC"},
//    // 110000 x 110000
//    YMIR_ASSERT3(m.aaProbability("HI", 4, 4, 48, 48, 48), m_["CA"] * 2 + m_["TA"] * 2)
//    YMIR_ASSERT3(m.aaProbability("HI", 4, 5, 48, 48, 48), comp_prob("CAT", m_) * 2
//                                                          + comp_prob("TAT", m_) * 2)
//    YMIR_ASSERT3(m.aaProbability("HI", 4, 6, 48, 48, 48), comp_prob("CATT", m_)
//                                                          + comp_prob("CATC", m_)
//                                                          + comp_prob("TATC", m_)
//                                                          + comp_prob("TATT", m_))
//    // 100000 x 101000
//    YMIR_ASSERT3(m.aaProbability("HI", 4, 4, 40, 40, 32), 2 * m_["TA"])
//    YMIR_ASSERT3(m.aaProbability("HI", 4, 5, 40, 40, 32), 2 * comp_prob("TAT", m_))
//    YMIR_ASSERT3(m.aaProbability("HI", 4, 6, 40, 40, 32), comp_prob("TATT", m_) + comp_prob("TATA", m_))
//
//
//    //
//    // neighbour aminoacids
//    //
//    // {'A', "GCT"}, {'A', "GCC"}, {'A', "GCA"}, {'A', "GCG"}
//    // {'F', "TTT"}, {'F', "TTC"}
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 5, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 6, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("AF", 2, 4, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("AF", 2, 5, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("AF", 2, 6, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("AF", 3, 4, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("AF", 3, 5, 0, 0), 0)
//    YMIR_ASSERT3(m.aaProbability("AF", 3, 6, 0, 0), 0)
//
//    // 100000 and 100000
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 32, 32), .25 * comp_prob("GCTT", m_))
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 33, 32), .25 * comp_prob("GCTT", m_))
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 32, 33), .25 * comp_prob("GCTT", m_))
//
//    // 010000 and 010000
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 16, 16), .25 * comp_prob("GCCT", m_))
//
//    // 100000 and 010000
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 32, 16), .25 * comp_prob("GCTT", m_))
//
//    // 110000 and 100000
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 48, 32), .25 * (comp_prob("GCTT", m_) + comp_prob("GCCT", m_)))
//
//    // 100000 and 110000
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 32, 48), .25 * (comp_prob("GCTT", m_) * 2))
//
//    // 111100 and 110000 (full)
//    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 60, 48), .25 * 2 * (comp_prob("GCTT", m_)
//                                                                 + comp_prob("GCCT", m_)
//                                                                 + comp_prob("GCAT", m_)
//                                                                 + comp_prob("GCGT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 1, 5, 60, 48), .25 * 2 * (comp_prob("GCTTT", m_)
                                                                 + comp_prob("GCCTT", m_)
                                                                 + comp_prob("GCATT", m_)
                                                                 + comp_prob("GCGTT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 1, 6, 60, 48), .25 * (comp_prob("GCTTTT", m_)
                                                             + comp_prob("GCTTTC", m_)
                                                             + comp_prob("GCCTTT", m_)
                                                             + comp_prob("GCCTTC", m_)
                                                             + comp_prob("GCATTT", m_)
                                                             + comp_prob("GCATTC", m_)
                                                             + comp_prob("GCGTTT", m_)
                                                             + comp_prob("GCGTTC", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 2, 4, 60, 48), 2 * (comp_prob("GCTT", m_)
                                                           + comp_prob("GCCT", m_)
                                                           + comp_prob("GCAT", m_)
                                                           + comp_prob("GCGT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 2, 5, 60, 48), 2 * (comp_prob("GCTTT", m_)
                                                           + comp_prob("GCCTT", m_)
                                                           + comp_prob("GCATT", m_)
                                                           + comp_prob("GCGTT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 2, 6, 60, 48), (comp_prob("GCTTTT", m_)
                                                       + comp_prob("GCTTTC", m_)
                                                       + comp_prob("GCCTTT", m_)
                                                       + comp_prob("GCCTTC", m_)
                                                       + comp_prob("GCATTT", m_)
                                                       + comp_prob("GCATTC", m_)
                                                       + comp_prob("GCGTTT", m_)
                                                       + comp_prob("GCGTTC", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 3, 4, 60, 48), 2 * (comp_prob("CTT", m_)
                                                           + comp_prob("CCT", m_)
                                                           + comp_prob("CAT", m_)
                                                           + comp_prob("CGT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 3, 5, 60, 48), 2 * (comp_prob("CTTT", m_)
                                                           + comp_prob("CCTT", m_)
                                                           + comp_prob("CATT", m_)
                                                           + comp_prob("CGTT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 3, 6, 60, 48), (comp_prob("CTTTT", m_)
                                                       + comp_prob("CTTTC", m_)
                                                       + comp_prob("CCTTT", m_)
                                                       + comp_prob("CCTTC", m_)
                                                       + comp_prob("CATTT", m_)
                                                       + comp_prob("CATTC", m_)
                                                       + comp_prob("CGTTT", m_)
                                                       + comp_prob("CGTTC", m_)))


    // distant aminoacids
    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_markovchain_aa_di_rev)

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_markovchain_nuc_mono_err)

    vector<prob_t> probs = {.1, .2, .3, .4};
    MonoNucInsertionModel m(probs.begin(), .5);
    YMIR_ASSERT2(m.err_prob(), .5)

    string s = "ACGT";

    // (.1 + .5 * (.2 + .3 + .4))
    // * (.2 + .5 * (.1 + .3 + .4))
    // * (.3 + .5 * (.2 + .1 + .4))
    // * (.4 + .5 * (.2 + .3 + .1))
    YMIR_ASSERT(abs(m.nucProbability(s.begin(), 4, NULL_CHAR, true) - 0.15015) < 1e-17)
    YMIR_ASSERT(abs(m.nucProbability(s, NULL_CHAR, true) - 0.15015) < 1e-17)

    // (.1 + .5 * (.2 + .3 + .4))
    // * (.2 + .5 * (.1 + .3 + .4))
    // * (.3 + .5 * (.2 + .1 + .4))
    // * (.4 + .5 * (.2 + .3 + .1))
    YMIR_ASSERT(abs(m.nucProbability(s.begin(), 4, 'A', true) - 0.15015) < 1e-17)
    YMIR_ASSERT(abs(m.nucProbability(s, 'A', true) - 0.15015) < 1e-17)

    // (.1 + .5 * (.2 + .3 + .4))
    // * (.2 + .5 * (.1 + .3 + .4))
    // * (.3 + .5 * (.2 + .1 + .4))
    YMIR_ASSERT(abs(m.nucProbability(s.substr(0, 3), 'A', true) - 0.2145) < 1e-17)
    YMIR_ASSERT(abs(m.nucProbability(s.substr(0, 3), NULL_CHAR, true) - 0.2145) < 1e-17)

    probs = {1, 0, 0, 0};
    m = MonoNucInsertionModel(probs.begin());
    std::default_random_engine rg;
    YMIR_ASSERT2(m.generate(5, rg), "AAAAA");

    probs = {0, 0, 0, 1};
    m = MonoNucInsertionModel(probs.begin());
    YMIR_ASSERT2(m.generate(1, rg), "T");

    probs = {0, 1, 0, 0};
    m = MonoNucInsertionModel(probs.begin());
    YMIR_ASSERT2(m.generate(3, rg), "CCC");

YMIR_TEST_END


YMIR_TEST_START(test_markovchain_nuc_di_err)

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

    DiNucInsertionModel mc(mat, .5);
    YMIR_ASSERT2(mc.err_prob(), .5)

    string s = "ACGT";

    // err:    .25 * (.25 + .75 * .7) * (.25 + .75 * .3) * (.25 + .75 * .2)
    // no err: .25 * .7 * .3 * .2 = .0105
    YMIR_ASSERT(mc.nucProbability("", NULL_CHAR, true) == 1);
    YMIR_ASSERT(abs(mc.nucProbability(s, NULL_CHAR, true) - 0.0368125) < 1e-15);
    YMIR_ASSERT(abs(mc.nucProbability(s.begin(), 4, NULL_CHAR, true) - 0.0368125) < 1e-15);

    // err:    (.25 + .75 * .4) * (.25 + .75 * .7) * (.25 + .75 * .3) * (.25 + .75 * .2)
    // no err: .4 * .7 * .3 * .2 = .0126
    YMIR_ASSERT(abs(mc.nucProbability(s.begin(), 4, 'C', true) - 0.0809875) < 1e-17);

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

    // Tests for markov chain.
//    YMIR_TEST(test_markovchain_nuc_mono())
//    YMIR_TEST(test_markovchain_nuc_di())
//    YMIR_TEST(test_markovchain_aa_mono())
    YMIR_TEST(test_markovchain_aa_di())
//    YMIR_TEST(test_markovchain_aa_di_rev())
//    YMIR_TEST(test_markovchain_nuc_mono_err())
//    YMIR_TEST(test_markovchain_nuc_di_err())

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