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

    auto aapr = [](const std::string &str, const  std::map<std::string, prob_t> &m_) {
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
    YMIR_ASSERT3(m.aaProbability("M", 1, 1, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("M", 1, 2, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("M", 1, 3, 0, 0), 0)

    YMIR_ASSERT3(m.aaProbability("M", 0, 1, 63, 63), .25)
    YMIR_ASSERT3(m.aaProbability("M", 1, 1, 63, 63), .25)
    YMIR_ASSERT3(m.aaProbability("M", 1, 2, 63, 63), .25 * .4)
    YMIR_ASSERT3(m.aaProbability("M", 1, 3, 63, 63), .25 * .4 * .3)
    YMIR_ASSERT3(m.aaProbability("M", 2, 2, 63, 63), .4)
    YMIR_ASSERT3(m.aaProbability("M", 2, 3, 63, 63), .4 * .3)
    YMIR_ASSERT3(m.aaProbability("M", 3, 3, 63, 63), .3)

    // {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 0, 0), 0)

    // 111000
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 56, 56), m_["GA"] + m_["GA"] + m_["GA"])
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 56, 56), aapr("GAT", m_) * 3)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 56, 56), aapr("GATT", m_) + aapr("GATC", m_) + aapr("GATA", m_))
    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 56, 56), m_["AT"] + m_["AT"] + m_["AT"])
    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 56, 56), aapr("ATT", m_) + aapr("ATC", m_) + aapr("ATA", m_))
    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 56, 56), m_["TT"] + m_["TC"] + m_["TA"])

    // 101000
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 40, 40), m_["GA"] + m_["GA"])
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 40, 40), aapr("GAT", m_) * 2)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 40, 40), aapr("GATT", m_) + aapr("GATA", m_))
    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 40, 40), m_["AT"] + m_["AT"])
    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 40, 40), aapr("ATT", m_) + aapr("ATA", m_))
    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 40, 40), m_["TT"] + m_["TA"])

    // 110000 and 011000 -> 010000
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 48, 24), m_["GA"])
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 48, 24), aapr("GAT", m_))
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 48, 24), aapr("GATC", m_))
    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 48, 24), m_["AT"])
    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 48, 24), aapr("ATC", m_))
    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 48, 24), m_["TC"])

    // 100000 and 101000 -> 100000
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 32, 40), m_["GA"])
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 48, 40), aapr("GAT", m_))
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 48, 40), aapr("GATT", m_))

    // prev: {'H', "CAT"}, {'H', "CAC"},
    // 110000 x 110000
    YMIR_ASSERT3(m.aaProbability("HI", 4, 4, 48, 48, 48), m_["CA"] * 2 + m_["TA"] * 2)
    YMIR_ASSERT3(m.aaProbability("HI", 4, 5, 48, 48, 48), aapr("CAT", m_) * 2
                                                          + aapr("TAT", m_) * 2)
    YMIR_ASSERT3(m.aaProbability("HI", 4, 6, 48, 48, 48), aapr("CATT", m_)
                                                          + aapr("CATC", m_)
                                                          + aapr("TATC", m_)
                                                          + aapr("TATT", m_))
    // 100000 x 101000
    YMIR_ASSERT3(m.aaProbability("HI", 4, 4, 40, 40, 32), 2 * m_["TA"])
    YMIR_ASSERT3(m.aaProbability("HI", 4, 5, 40, 40, 32), 2 * aapr("TAT", m_))
    YMIR_ASSERT3(m.aaProbability("HI", 4, 6, 40, 40, 32), aapr("TATT", m_) + aapr("TATA", m_))


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

    // 100000 and 100000
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 32, 32), .25 * aapr("GCTT", m_))
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 33, 32), .25 * aapr("GCTT", m_))
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 32, 33), .25 * aapr("GCTT", m_))

    // 010000 and 010000
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 16, 16), .25 * aapr("GCCT", m_))

    // 100000 and 010000
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 32, 16), .25 * aapr("GCTT", m_))

    // 110000 and 100000
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 48, 32), .25 * (aapr("GCTT", m_) + aapr("GCCT", m_)))

    // 100000 and 110000
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 32, 48), .25 * (aapr("GCTT", m_) * 2))

    // 111100 and 110000 (full)
    YMIR_ASSERT3(m.aaProbability("AF", 1, 4, 60, 48), .25 * 2 * (aapr("GCTT", m_)
                                                                 + aapr("GCCT", m_)
                                                                 + aapr("GCAT", m_)
                                                                 + aapr("GCGT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 1, 5, 60, 48), .25 * 2 * (aapr("GCTTT", m_)
                                                                 + aapr("GCCTT", m_)
                                                                 + aapr("GCATT", m_)
                                                                 + aapr("GCGTT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 1, 6, 60, 48), .25 * (aapr("GCTTTT", m_)
                                                             + aapr("GCTTTC", m_)
                                                             + aapr("GCCTTT", m_)
                                                             + aapr("GCCTTC", m_)
                                                             + aapr("GCATTT", m_)
                                                             + aapr("GCATTC", m_)
                                                             + aapr("GCGTTT", m_)
                                                             + aapr("GCGTTC", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 2, 4, 60, 48), 2 * (aapr("GCTT", m_)
                                                           + aapr("GCCT", m_)
                                                           + aapr("GCAT", m_)
                                                           + aapr("GCGT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 2, 5, 60, 48), 2 * (aapr("GCTTT", m_)
                                                           + aapr("GCCTT", m_)
                                                           + aapr("GCATT", m_)
                                                           + aapr("GCGTT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 2, 6, 60, 48), (aapr("GCTTTT", m_)
                                                       + aapr("GCTTTC", m_)
                                                       + aapr("GCCTTT", m_)
                                                       + aapr("GCCTTC", m_)
                                                       + aapr("GCATTT", m_)
                                                       + aapr("GCATTC", m_)
                                                       + aapr("GCGTTT", m_)
                                                       + aapr("GCGTTC", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 3, 4, 60, 48), 2 * (aapr("CTT", m_)
                                                           + aapr("CCT", m_)
                                                           + aapr("CAT", m_)
                                                           + aapr("CGT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 3, 5, 60, 48), 2 * (aapr("CTTT", m_)
                                                           + aapr("CCTT", m_)
                                                           + aapr("CATT", m_)
                                                           + aapr("CGTT", m_)))
    YMIR_ASSERT3(m.aaProbability("AF", 3, 6, 60, 48), (aapr("CTTTT", m_)
                                                       + aapr("CTTTC", m_)
                                                       + aapr("CCTTT", m_)
                                                       + aapr("CCTTC", m_)
                                                       + aapr("CATTT", m_)
                                                       + aapr("CATTC", m_)
                                                       + aapr("CGTTT", m_)
                                                       + aapr("CGTTC", m_)))


    // distant aminoacids
    // {'H', "CAT"}, {'H', "CAC"},
    // {'F', "TTT"}, {'F', "TTC"}
    // {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
    // 100000 x 110000 x 001000
    YMIR_ASSERT3(m.aaProbability("HFI", 1, 7, 32, 8), .25 * (  aapr("CATTTTA", m_)
                                                             + aapr("CATTTCA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 1, 8, 32, 8), .25 * (  aapr("CATTTTAT", m_)
                                                               + aapr("CATTTCAT", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 1, 9, 32, 8), .25 * (  aapr("CATTTTATA", m_)
                                                               + aapr("CATTTCATA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 2, 7, 32, 8), (  aapr("CATTTTA", m_)
                                                         + aapr("CATTTCA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 2, 8, 32, 8), (  aapr("CATTTTAT", m_)
                                                         + aapr("CATTTCAT", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 2, 9, 32, 8), (  aapr("CATTTTATA", m_)
                                                         + aapr("CATTTCATA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 3, 7, 32, 8),   (aapr("ATTTTA", m_)
                                                       + aapr("ATTTCA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 3, 8, 32, 8),   (aapr("ATTTTAT", m_)
                                                       + aapr("ATTTCAT", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 3, 9, 32, 8),   (aapr("ATTTTATA", m_)
                                                       + aapr("ATTTCATA", m_)))

    // 110000 x 110000 x 001000
    YMIR_ASSERT3(m.aaProbability("HFI", 1, 7, 48, 8), .25 * (  aapr("CATTTTA", m_)
                                                             + aapr("CACTTTA", m_)
                                                             + aapr("CATTTCA", m_)
                                                             + aapr("CACTTCA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 1, 8, 48, 8), .25 * (  aapr("CATTTTAT", m_)
                                                             + aapr("CACTTTAT", m_)
                                                             + aapr("CATTTCAT", m_)
                                                             + aapr("CACTTCAT", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 1, 9, 48, 8), .25 * (  aapr("CATTTTATA", m_)
                                                             + aapr("CACTTTATA", m_)
                                                             + aapr("CATTTCATA", m_)
                                                             + aapr("CACTTCATA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 2, 7, 48, 8), (  aapr("CATTTTA", m_)
                                                       + aapr("CACTTTA", m_)
                                                       + aapr("CATTTCA", m_)
                                                       + aapr("CACTTCA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 2, 8, 48, 8), (  aapr("CATTTTAT", m_)
                                                       + aapr("CACTTTAT", m_)
                                                       + aapr("CATTTCAT", m_)
                                                       + aapr("CACTTCAT", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 2, 9, 48, 8), (  aapr("CATTTTATA", m_)
                                                       + aapr("CACTTTATA", m_)
                                                       + aapr("CATTTCATA", m_)
                                                       + aapr("CACTTCATA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 3, 7, 48, 8), (  aapr("ATTTTA", m_)
                                                       + aapr("ACTTTA", m_)
                                                       + aapr("ATTTCA", m_)
                                                       + aapr("ACTTCA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 3, 8, 48, 8), (  aapr("ATTTTAT", m_)
                                                       + aapr("ACTTTAT", m_)
                                                       + aapr("ATTTCAT", m_)
                                                       + aapr("ACTTCAT", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 3, 9, 48, 8), (  aapr("ATTTTATA", m_)
                                                       + aapr("ACTTTATA", m_)
                                                       + aapr("ATTTCATA", m_)
                                                       + aapr("ACTTCATA", m_)))

    // {'H', "CAT"}, {'H', "CAC"},
    // {'F', "TTT"}, {'F', "TTC"}
    // {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
    // 100000 x 110000 x 011000
    YMIR_ASSERT3(m.aaProbability("HFI", 1, 7, 32, 24), 2 * .25 * (  aapr("CATTTTA", m_)
                                                                  + aapr("CATTTCA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 1, 8, 32, 24), 2 * .25 * (  aapr("CATTTTAT", m_)
                                                                  + aapr("CATTTCAT", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 1, 9, 32, 24), .25 * (  aapr("CATTTTATA", m_)
                                                              + aapr("CATTTCATA", m_)
                                                              + aapr("CATTTTATC", m_)
                                                              + aapr("CATTTCATC", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 2, 7, 32, 24), 2 * (aapr("CATTTTA", m_)
                                                          + aapr("CATTTCA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 2, 8, 32, 24), 2 * (aapr("CATTTTAT", m_)
                                                          + aapr("CATTTCAT", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 2, 9, 32, 24), (  aapr("CATTTTATA", m_)
                                                        + aapr("CATTTCATA", m_)
                                                        + aapr("CATTTTATC", m_)
                                                        + aapr("CATTTCATC", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 3, 7, 32, 24), 2 * (aapr("ATTTTA", m_)
                                                          + aapr("ATTTCA", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 3, 8, 32, 24), 2 * (aapr("ATTTTAT", m_)
                                                          + aapr("ATTTCAT", m_)))
    YMIR_ASSERT3(m.aaProbability("HFI", 3, 9, 32, 24), ( aapr("ATTTTATA", m_)
                                                       + aapr("ATTTCATA", m_)
                                                       + aapr("ATTTTATC", m_)
                                                       + aapr("ATTTCATC", m_)))

YMIR_TEST_END


YMIR_TEST_START(test_markovchain_aa_di_rev)

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

    // note: backward probabilbities
    // prev -> next
    std::map<std::string, prob_t> m_;
    m_["AA"] = mat(0, 0);
    m_["AC"] = mat(1, 0);
    m_["AG"] = mat(2, 0);
    m_["AT"] = mat(3, 0);

    m_["CA"] = mat(0, 1);
    m_["CC"] = mat(1, 1);
    m_["CG"] = mat(2, 1);
    m_["CT"] = mat(3, 1);

    m_["GA"] = mat(0, 2);
    m_["GC"] = mat(1, 2);
    m_["GG"] = mat(2, 2);
    m_["GT"] = mat(3, 2);

    m_["TA"] = mat(0, 3);
    m_["TC"] = mat(1, 3);
    m_["TG"] = mat(2, 3);
    m_["TT"] = mat(3, 3);

    auto aapr = [](const std::string &str, const  std::map<std::string, prob_t> &m_) {
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
    YMIR_ASSERT3(m.aaProbabilityRev("M", 1, 1, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("M", 2, 1, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("M", 3, 1, 0, 0), 0)

    YMIR_ASSERT3(m.aaProbabilityRev("M", 4, 3, 63, 63), .25)
    YMIR_ASSERT3(m.aaProbabilityRev("M", 3, 3, 63, 63), .25)
    YMIR_ASSERT3(m.aaProbabilityRev("M", 3, 2, 63, 63), .25 * m_["TG"])
    YMIR_ASSERT3(m.aaProbabilityRev("M", 3, 1, 63, 63), .25 * m_["TG"] * m_["AT"])
    YMIR_ASSERT3(m.aaProbabilityRev("M", 2, 2, 63, 63), m_["TG"])
    YMIR_ASSERT3(m.aaProbabilityRev("M", 2, 1, 63, 63), m_["TG"] * m_["AT"])
    YMIR_ASSERT3(m.aaProbabilityRev("M", 1, 1, 63, 63), m_["AT"])

    // {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 6, 6, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 6, 5, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 6, 4, 0, 0), 0)

    // 111000
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 3, 56, 56), aapr("TA", m_) + aapr("CA", m_) + aapr("AA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 2, 56, 56), aapr("TTA", m_) + aapr("TCA", m_) + aapr("TAA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 1, 56, 56), aapr("ATTA", m_) + aapr("ATCA", m_) + aapr("ATAA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 2, 2, 56, 56), aapr("TT", m_) + aapr("TC", m_) + aapr("TA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 2, 1, 56, 56), aapr("ATT", m_) + aapr("ATC", m_) + aapr("ATA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 1, 1, 56, 56), 3 * aapr("AT", m_))

    // 101000
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 3, 40, 40), aapr("AT", m_) + aapr("AA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 2, 40, 40), aapr("TTA", m_) + aapr("TAA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 1, 40, 40), aapr("ATTA", m_) + aapr("ATAA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 2, 2, 40, 40), aapr("TT", m_) + aapr("TA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 2, 1, 40, 40), aapr("ATT", m_) + aapr("ATA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 1, 1, 40, 40), 2 * aapr("TA", m_))

    // 110000 and 011000 -> 010000
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 3, 48, 24), aapr("CA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 2, 48, 24), aapr("TCA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 1, 48, 24), aapr("ATCA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 2, 2, 48, 24), aapr("TC", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 2, 1, 48, 24), aapr("ATC", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 1, 1, 48, 24), aapr("AT", m_))

    // 100000 and 101000 -> 100000
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 3, 32, 40), aapr("TA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 2, 32, 40), aapr("TTA", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IM", 3, 1, 32, 40), aapr("ATTA", m_))

    // prev: {'H', "CAT"}, {'H', "CAC"},
    // next: {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
    // 110000 x 110000
    YMIR_ASSERT3(m.aaProbabilityRev("IH", 3, 3, 48, 48, 48), 2 * (aapr("TC", m_) + aapr("CC", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IH", 3, 2, 48, 48, 48), 2 * (aapr("TTC", m_) + aapr("TCT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IH", 3, 1, 48, 48, 48), 2 * aapr("ATTC", m_)
                                                             + 2 * aapr("ATCC", m_))
    // (prev, i.e., H) 100000 x 101000 (next, i.e., I)
    YMIR_ASSERT3(m.aaProbabilityRev("IH", 3, 3, 40, 40, 32), aapr("TC", m_) + aapr("AC", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IH", 3, 2, 40, 40, 32), aapr("TTC", m_) + aapr("TAC", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("IH", 3, 1, 40, 40, 32), aapr("ATTC", m_) + aapr("ATAC", m_))


    //
    // neighbour aminoacids
    //
    // {'F', "TTT"}, {'F', "TTC"}
    // {'A', "GCT"}, {'A', "GCC"}, {'A', "GCA"}, {'A', "GCG"}
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 3, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 2, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 1, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 5, 3, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 5, 2, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 5, 3, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 4, 1, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 4, 2, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 4, 1, 0, 0), 0)

    // 100000 and 100000
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 3, 32, 32), .25 * aapr("TGCT", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 3, 33, 32), .25 * aapr("TGCT", m_))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 3, 32, 33), .25 * aapr("TGCT", m_))

    // 010000 and 010000
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 3, 16, 16), .25 * aapr("CGCC", m_))

    // 010000 and 100000
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 3, 32, 16), .25 * aapr("CGCC", m_))

    // 100000 and 110000
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 3, 48, 32), .25 * (aapr("TGCT", m_) + aapr("TGCC", m_)))

    // 110000 and 100000
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 3, 32, 48), .25 * (aapr("TGCT", m_) + aapr("CGCT", m_)))

    // A:111100 and F:110000 (full)
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 3, 60, 48), .25 * (aapr("TGCT", m_) + aapr("CGCT", m_)
                                                                 +  aapr("TGCC", m_) + aapr("CGCC", m_)
                                                                 +  aapr("TGCA", m_) + aapr("CGCA", m_)
                                                                 +  aapr("TGCG", m_) + aapr("CGCG", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 2, 60, 48), .25 * (aapr("TTGCT", m_)
                                                             +  aapr("TTGCC", m_)
                                                             +  aapr("TTGCA", m_)
                                                             +  aapr("TTGCG", m_)
                                                             +  aapr("TCGCT", m_)
                                                             +  aapr("TCGCC", m_)
                                                             +  aapr("TCGCA", m_)
                                                             +  aapr("TCGCG", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 6, 1, 60, 48), .25 * (aapr("TTTGCT", m_)
                                                             +  aapr("TTTGCC", m_)
                                                             +  aapr("TTTGCA", m_)
                                                             +  aapr("TTTGCG", m_)
                                                             +  aapr("TTCGCT", m_)
                                                             +  aapr("TTCGCC", m_)
                                                             +  aapr("TTCGCA", m_)
                                                             +  aapr("TTCGCG", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 5, 3, 60, 48), (aapr("TGCT", m_) + aapr("CGCT", m_)
                                                          +  aapr("TGCC", m_) + aapr("CGCC", m_)
                                                          +  aapr("TGCA", m_) + aapr("CGCA", m_)
                                                          +  aapr("TGCG", m_) + aapr("CGCG", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 5, 2, 60, 48),    (aapr("TTGCT", m_)
                                                          +  aapr("TTGCC", m_)
                                                          +  aapr("TTGCA", m_)
                                                          +  aapr("TTGCG", m_)
                                                          +  aapr("TCGCT", m_)
                                                          +  aapr("TCGCC", m_)
                                                          +  aapr("TCGCA", m_)
                                                          +  aapr("TCGCG", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 5, 1, 60, 48),    (aapr("TTTGCT", m_)
                                                          +  aapr("TTTGCC", m_)
                                                          +  aapr("TTTGCA", m_)
                                                          +  aapr("TTTGCG", m_)
                                                          +  aapr("TTCGCT", m_)
                                                          +  aapr("TTCGCC", m_)
                                                          +  aapr("TTCGCA", m_)
                                                          +  aapr("TTCGCG", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 4, 3, 60, 48), (aapr("TGC", m_) + aapr("CGC", m_)
                                                          +  aapr("TGC", m_) + aapr("CGC", m_)
                                                          +  aapr("TGC", m_) + aapr("CGC", m_)
                                                          +  aapr("TGC", m_) + aapr("CGC", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 4, 2, 60, 48),    (aapr("TTGC", m_)
                                                          +  aapr("TTGC", m_)
                                                          +  aapr("TTGC", m_)
                                                          +  aapr("TTGC", m_)
                                                          +  aapr("TCGC", m_)
                                                          +  aapr("TCGC", m_)
                                                          +  aapr("TCGC", m_)
                                                          +  aapr("TCGC", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("FA", 4, 1, 60, 48),   (aapr("TTTGC", m_)
                                                         +  aapr("TTTGC", m_)
                                                         +  aapr("TTTGC", m_)
                                                         +  aapr("TTTGC", m_)
                                                         +  aapr("TTCGC", m_)
                                                         +  aapr("TTCGC", m_)
                                                         +  aapr("TTCGC", m_)
                                                         +  aapr("TTCGC", m_)))


    // distant aminoacids
    // {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
    // {'F', "TTT"}, {'F', "TTC"}
    // {'H', "CAT"}, {'H', "CAC"},
    // I:001000 x F:110000 x H:100000
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 4, 32, 63), .25 * (aapr("TTTCAT", m_)
                                                               + aapr("TTCCAT", m_)))

    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 3, 32, 8), .25 * (aapr("ATTTCAT", m_)
                                                              + aapr("ATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 2, 32, 8), .25 * (aapr("TATTTCAT", m_)
                                                              + aapr("TATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 1, 32, 8), .25 * (aapr("ATATTTCAT", m_)
                                                              + aapr("ATATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 8, 3, 32, 8), (aapr("ATTTCAT", m_)
                                                          + aapr("ATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 8, 2, 32, 8), (aapr("TATTTCAT", m_)
                                                          + aapr("TATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 8, 1, 32, 8), (aapr("ATATTTCAT", m_)
                                                          + aapr("ATATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 7, 3, 32, 8), (aapr("ATTTCA", m_)
                                                          + aapr("ATTCCA", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 7, 2, 32, 8), (aapr("TATTTCA", m_)
                                                          + aapr("TATTCCA", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 7, 1, 32, 8), (aapr("ATATTTCA", m_)
                                                          + aapr("ATATTCCA", m_)))

    //  001000 x 110000 x 110000
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 3, 48, 8), .25 * ( aapr("ATTTCAT", m_)
                                                               + aapr("ATTCCAT", m_)
                                                               + aapr("ATTTCAC", m_)
                                                               + aapr("ATTCCAC", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 2, 48, 8), .25 * ( aapr("TATTTCAT", m_)
                                                               + aapr("TATTCCAT", m_)
                                                               + aapr("TATTTCAC", m_)
                                                               + aapr("TATTCCAC", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 1, 48, 8), .25 * ( aapr("ATATTTCAT", m_)
                                                               + aapr("ATATTCCAT", m_)
                                                               + aapr("ATATTTCAC", m_)
                                                               + aapr("ATATTCCAC", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 8, 3, 48, 8), ( aapr("ATTTCAT", m_)
                                                         + aapr("ATTCCAT", m_)
                                                         + aapr("ATTTCAC", m_)
                                                         + aapr("ATTCCAC", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 8, 2, 48, 8), ( aapr("TATTTCAT", m_)
                                                         + aapr("TATTCCAT", m_)
                                                         + aapr("TATTTCAC", m_)
                                                         + aapr("TATTCCAC", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 8, 1, 48, 8), ( aapr("ATATTTCAT", m_)
                                                         + aapr("ATATTCCAT", m_)
                                                         + aapr("ATATTTCAC", m_)
                                                         + aapr("ATATTCCAC", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 7, 3, 48, 8), ( aapr("ATTTCA", m_)
                                                         + aapr("ATTCCA", m_)
                                                         + aapr("ATTTCA", m_)
                                                         + aapr("ATTCCA", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 7, 2, 48, 8), ( aapr("TATTTCA", m_)
                                                         + aapr("TATTCCA", m_)
                                                         + aapr("TATTTCA", m_)
                                                         + aapr("TATTCCA", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 7, 1, 48, 8), ( aapr("ATATTTCA", m_)
                                                         + aapr("ATATTCCA", m_)
                                                         + aapr("ATATTTCA", m_)
                                                         + aapr("ATATTCCA", m_)))

    // {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
    // {'F', "TTT"}, {'F', "TTC"}
    // {'H', "CAT"}, {'H', "CAC"},
    // 011000 x 110000 x 100000
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 3, 32, 24), .25 * ( aapr("CTTTCAT", m_)
                                                                + aapr("CTTCCAT", m_)
                                                                + aapr("ATTTCAT", m_)
                                                                + aapr("ATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 2, 32, 24), .25 * ( aapr("TCTTTCAT", m_)
                                                                + aapr("TCTTCCAT", m_)
                                                                + aapr("TATTTCAT", m_)
                                                                + aapr("TATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 9, 1, 32, 24), .25 * ( aapr("ATCTTTCAT", m_)
                                                                + aapr("ATCTTCCAT", m_)
                                                                + aapr("ATATTTCAT", m_)
                                                                + aapr("ATATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 8, 3, 32, 24), ( aapr("CTTTCAT", m_)
                                                          + aapr("CTTCCAT", m_)
                                                          + aapr("ATTTCAT", m_)
                                                          + aapr("ATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 8, 2, 32, 24), ( aapr("TCTTTCAT", m_)
                                                          + aapr("TCTTCCAT", m_)
                                                          + aapr("TATTTCAT", m_)
                                                          + aapr("TATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 8, 1, 32, 24), ( aapr("ATCTTTCAT", m_)
                                                          + aapr("ATCTTCCAT", m_)
                                                          + aapr("ATATTTCAT", m_)
                                                          + aapr("ATATTCCAT", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 7, 3, 32, 24), ( aapr("CTTTCA", m_)
                                                          + aapr("CTTCCA", m_)
                                                          + aapr("ATTTCA", m_)
                                                          + aapr("ATTCCA", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 7, 2, 32, 24), ( aapr("TCTTTCA", m_)
                                                          + aapr("TCTTCCA", m_)
                                                          + aapr("TATTTCA", m_)
                                                          + aapr("TATTCCA", m_)))
    YMIR_ASSERT3(m.aaProbabilityRev("IFH", 7, 1, 32, 24), ( aapr("ATCTTTCA", m_)
                                                          + aapr("ATCTTCCA", m_)
                                                          + aapr("ATATTTCA", m_)
                                                          + aapr("ATATTCCA", m_)))

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
    YMIR_TEST(test_markovchain_nuc_mono())
    YMIR_TEST(test_markovchain_nuc_di())
    YMIR_TEST(test_markovchain_aa_mono())
    YMIR_TEST(test_markovchain_aa_di())
    YMIR_TEST(test_markovchain_aa_di_rev())
    YMIR_TEST(test_markovchain_nuc_mono_err())
    YMIR_TEST(test_markovchain_nuc_di_err())

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