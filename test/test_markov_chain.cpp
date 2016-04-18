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

    // exact same aminoacid
    // {'M', "ATG"}
    YMIR_ASSERT3(m.aaProbability("M", 1, 1, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("M", 1, 2, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("M", 1, 3, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("M", 1, 1, 63, 63), .1)
    YMIR_ASSERT3(m.aaProbability("M", 1, 2, 63, 63), .1 * .4)
    YMIR_ASSERT3(m.aaProbability("M", 1, 3, 63, 63), .1 * .4 * .3)
    YMIR_ASSERT3(m.aaProbability("M", 2, 3, 63, 63), .4 * .3)
    YMIR_ASSERT3(m.aaProbability("M", 3, 3, 63, 63), .3)

    // {'I', "ATT"}, {'I', "ATC"}, {'I', "ATA"}
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 0, 0), 0)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 0, 0), 0)

    // 111000
    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 56, 56), 3 * .1)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 56, 56), 3 * .1 * 3 * .4)
    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 56, 56), (3 * .1) * (3 * .4) * (.1 + .2 + .4))
    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 56, 56), (3 * .4))
    YMIR_ASSERT3(m.aaProbability("MI", 5, 6, 56, 56), (3 * .4) * (.1 + .2 + .4))
    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 56, 56), (.1 + .2 + .4))

    // 101000
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 40, 40), .1)
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 40, 40), .1 * .4)
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 40, 40), .1 * .4 * .3)
//    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 40, 40), .4 * .3)
//    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 40, 40), .3)
//
//    // 110000 and 011000
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 4, 48, 24), .1)
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 5, 48, 24), .1 * .4)
//    YMIR_ASSERT3(m.aaProbability("MI", 4, 6, 48, 24), .1 * .4 * .3)
//    YMIR_ASSERT3(m.aaProbability("MI", 5, 5, 48, 24), .4 * .3)
//    YMIR_ASSERT3(m.aaProbability("MI", 6, 6, 48, 24), .3)


    // neighbour aminoacids

    // distant aminoacids


    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_markovchain_aa_di)

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


    // neighbour aminoacids

    // exact same aminoacid

    // distant aminoacids


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
    YMIR_TEST(test_markovchain_nuc_mono())
    YMIR_TEST(test_markovchain_nuc_di())
    YMIR_TEST(test_markovchain_aa_mono())
    YMIR_TEST(test_markovchain_aa_di())
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