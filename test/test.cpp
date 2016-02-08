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

#define YMIR_TEST(res) { all_tests += 1; if (res.size() == 0) {tests_passed += 1;} \
                            else { \
                            failed_test_info.push_back(TestInfo(#res, res)); \
                            } \
                            }

#define YMIR_TEST_START(funname) vector<string> funname() { std::cout << " --- " << (#funname) << std::endl; vector<string> _failed_cases;
#define YMIR_ASSERT(expr) { if (!(expr)) { _failed_cases.push_back(#expr); } };
#define YMIR_ASSERT2(expr1, expr2) { if ((expr1) != (expr2)) { std::stringstream ss; ss << #expr1 << " == " << #expr2 << "  (result: " << (expr1) << ", need: " << (expr2) << ")";_failed_cases.push_back(ss.str()); } };
#define YMIR_TEST_END std::cout << std::endl; return _failed_cases; }

#define YMIR_TEST_PRECISION 1e-14

#include <iostream>
#include <list>
#include <sstream>

#include "Inference"


using namespace std;
using namespace ymir;


std::string TEST_DATA_FOLDER;


ModelParameterVector make_test_events_vj() {
    vector<prob_t> v1;  // param vec
    vector<event_ind_t> v2;  // lens vec
    vector<event_ind_t> v3;  // event classes
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
    v1.push_back(.1); v1.push_back(.2); v1.push_back(.3); v1.push_back(.4);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(8);

    return ModelParameterVector(VJ_RECOMB, v1, v2, v3, v4);
}


ModelParameterVector make_test_events_vj2() {
    vector<prob_t> v1;  // param vec
    vector<event_ind_t> v2;  // lens vec
    vector<event_ind_t> v3;  // event classes
    vector<seq_len_t> v4;  // event family col numbers

    // V-J
    v1.push_back(.049); v1.push_back(.026); v1.push_back(.035); // J1
    v1.push_back(.046); v1.push_back(.054); v1.push_back(.065); // J2
    v1.push_back(.075); v1.push_back(.084); v1.push_back(.566); // J3
    v2.push_back(9);

    v3.push_back(0);
    v4.push_back(3);

    // V del
    v1.push_back(.5); v1.push_back(.4); v1.push_back(.06); v1.push_back(.01); v1.push_back(.03);
    v2.push_back(5);

    v1.push_back(.4); v1.push_back(.11); v1.push_back(.1); v1.push_back(.39);
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
    v1.push_back(.1); v1.push_back(.05); v1.push_back(.2); v1.push_back(.15); v1.push_back(.24); v1.push_back(.01); v1.push_back(.25);
    v2.push_back(7);
    v4.push_back(0);

    v3.push_back(7);

    // VJ ins nuc
    v1.push_back(.4); v1.push_back(.3); v1.push_back(.2); v1.push_back(.1);
    v2.push_back(4);
    v4.push_back(0);

    v3.push_back(8);

    return ModelParameterVector(VJ_RECOMB, v1, v2, v3, v4);
}


ModelParameterVector make_test_events_vdj() {
    vector<prob_t> v1;  // param vec
    vector<event_ind_t> v2;  // lens vec
    vector<event_ind_t> v3;  // event classes
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
    // max 1 D3' del
    // max 1 D5' del
    v1.push_back(.17); // first row
    v1.push_back(.27);
    v1.push_back(.37); // second row
    v1.push_back(.19);

    v2.push_back(4);
    v4.push_back(2);

    // D2 dels
    // max 2 D3' del
    // max 1 D5' del
    // 3 rows 2 columns
    v1.push_back(.11); v1.push_back(.12);
    v1.push_back(.13); v1.push_back(.14);
    v1.push_back(.15); v1.push_back(.35);

    v2.push_back(6);
    v4.push_back(2);

    // D3 dels
    // max 1 D3' del
    // max 1 D5' del
    v1.push_back(.11); // first row
    v1.push_back(.22);
    v1.push_back(.33); // second row
    v1.push_back(.34);

    v2.push_back(4);
    v4.push_back(2);

    v3.push_back(8);

    // VD ins len
    v1.push_back(.05); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2); v1.push_back(.25); v1.push_back(.24); v1.push_back(.01);
    v2.push_back(7);

    v3.push_back(11);
    v4.push_back(0);

    // DJ ins len
    v1.push_back(.1); v1.push_back(.24); v1.push_back(.25); v1.push_back(.05); v1.push_back(.01); v1.push_back(.15); v1.push_back(.2);
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
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.15); v1.push_back(.5);
    v2.push_back(4);

    v3.push_back(15);
    v4.push_back(0);

    // prev T
    v1.push_back(.15); v1.push_back(.1); v1.push_back(.25); v1.push_back(.5);
    v2.push_back(4);

    v3.push_back(16);
    v4.push_back(0);

    // DJ ins nuc
    // prev A
    v1.push_back(.25); v1.push_back(.1); v1.push_back(.15); v1.push_back(.5);
    v2.push_back(4);

    v3.push_back(17);
    v4.push_back(0);

    // prev C
    v1.push_back(.15); v1.push_back(.1); v1.push_back(.25); v1.push_back(.5);
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

    vector<seq_len_t> v5; // min D genes len == 3
    v5.push_back(3);
    v5.push_back(3);
    v5.push_back(3);

    return ModelParameterVector(VDJ_RECOMB, v1, v2, v3, v4, vector<prob_t>(), true, v5);
}


ModelParameterVector make_test_events_vdj2() {
    vector<prob_t> v1;  // param vec
    vector<event_ind_t> v2;  // lens vec
    vector<event_ind_t> v3;  // event classes
    vector<seq_len_t> v4;  // event family col numbers

    // V
    v1.push_back(.5); v1.push_back(.15); v1.push_back(.35);
    v2.push_back(3);

    v3.push_back(0);
    v4.push_back(0);

    // J-D
    // J-D (3 Js - 3 Ds)
    v1.push_back(.02); v1.push_back(.03); v1.push_back(.03);
    v1.push_back(.03); v1.push_back(.04); v1.push_back(.06);
    v1.push_back(.07); v1.push_back(.1); v1.push_back(.62);

    v2.push_back(9);

    v3.push_back(1);
    v4.push_back(3);

    // V del
    v1.push_back(.5); v1.push_back(.4); v1.push_back(.04); v1.push_back(.03); v1.push_back(.03);
    v2.push_back(5);

    v1.push_back(.35); v1.push_back(.05); v1.push_back(.25); v1.push_back(.35);
    v2.push_back(4);

    v1.push_back(.74); v1.push_back(.006); v1.push_back(.02); v1.push_back(.01);
    v1.push_back(.04); v1.push_back(.039); v1.push_back(.145);
    v2.push_back(7);

    v3.push_back(2);
    v4.push_back(0);
    v4.push_back(0);
    v4.push_back(0);

    // J del
    v1.push_back(.08); v1.push_back(.1); v1.push_back(.2);
    v1.push_back(.33); v1.push_back(.13); v1.push_back(.15);
    v1.push_back(.01);
    v2.push_back(7);

    v1.push_back(.125); v1.push_back(.175); v1.push_back(.3);
    v1.push_back(.19); v1.push_back(.21);
    v2.push_back(5);

    v1.push_back(.11); v1.push_back(.1); v1.push_back(.02); v1.push_back(.02);
    v1.push_back(.02); v1.push_back(.03); v1.push_back(.7);
    v2.push_back(7);

    v3.push_back(5);
    v4.push_back(0);
    v4.push_back(0);
    v4.push_back(0);

    // D1 dels
    // max 1 D3' del
    // max 1 D5' del
    v1.push_back(.18); // first row
    v1.push_back(.28);
    v1.push_back(.36); // second row
    v1.push_back(.18);

    v2.push_back(4);
    v4.push_back(2);

    // D2 dels
    // max 2 D3' del
    // max 1 D5' del
    // 3 rows 2 columns
    v1.push_back(.12); v1.push_back(.11);
    v1.push_back(.14); v1.push_back(.13);
    v1.push_back(.16); v1.push_back(.34);

    v2.push_back(6);
    v4.push_back(2);

    // D3 dels
    // max 1 D3' del
    // max 1 D5' del
    v1.push_back(.12); // first row
    v1.push_back(.23);
    v1.push_back(.32); // second row
    v1.push_back(.33);

    v2.push_back(4);
    v4.push_back(2);

    v3.push_back(8);

    // VD ins len
    v1.push_back(.05); v1.push_back(.1); v1.push_back(.15); v1.push_back(.2); v1.push_back(.25); v1.push_back(.24); v1.push_back(.01);
    v2.push_back(7);

    v3.push_back(11);
    v4.push_back(0);

    // DJ ins len
    v1.push_back(.1); v1.push_back(.24); v1.push_back(.25); v1.push_back(.05); v1.push_back(.01); v1.push_back(.15); v1.push_back(.2);
    v2.push_back(7);

    v3.push_back(12);
    v4.push_back(0);

    // VD ins nuc
    // prev A
    v1.push_back(.04); v1.push_back(.09); v1.push_back(.03); v1.push_back(.84);
    v2.push_back(4);

    v3.push_back(13);
    v4.push_back(0);

    // prev C
    v1.push_back(.4); v1.push_back(.2); v1.push_back(.2); v1.push_back(.2);
    v2.push_back(4);

    v3.push_back(14);
    v4.push_back(0);

    // prev G
    v1.push_back(.25); v1.push_back(.15); v1.push_back(.1); v1.push_back(.5);
    v2.push_back(4);

    v3.push_back(15);
    v4.push_back(0);

    // prev T
    v1.push_back(.15); v1.push_back(.15); v1.push_back(.2); v1.push_back(.5);
    v2.push_back(4);

    v3.push_back(16);
    v4.push_back(0);

    // DJ ins nuc
    // prev A
    v1.push_back(.24); v1.push_back(.1); v1.push_back(.16); v1.push_back(.5);
    v2.push_back(4);

    v3.push_back(17);
    v4.push_back(0);

    // prev C
    v1.push_back(.16); v1.push_back(.1); v1.push_back(.24); v1.push_back(.5);
    v2.push_back(4);

    v3.push_back(18);
    v4.push_back(0);

    // prev G
    v1.push_back(.04); v1.push_back(.09); v1.push_back(.03); v1.push_back(.84);
    v2.push_back(4);

    v3.push_back(19);
    v4.push_back(0);

    // prev T
    v1.push_back(.4); v1.push_back(.11); v1.push_back(.29); v1.push_back(.2);
    v2.push_back(4);

    v3.push_back(20);
    v4.push_back(0);

    vector<seq_len_t> v5; // min D genes len == 3
    v5.push_back(3);
    v5.push_back(3);
    v5.push_back(3);

    return ModelParameterVector(VDJ_RECOMB, v1, v2, v3, v4, vector<prob_t>(), true, v5);
}


YMIR_TEST_START(test_basic)
YMIR_ASSERT(1 == 1)
YMIR_TEST_END


YMIR_TEST_START(test_model_param_vec_vj)
    vector<prob_t> v1;  // param vec
    vector<event_ind_t> v2;  // lens vec
    vector<event_ind_t> v3;  // event classes
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
    vector<event_ind_t> v2;  // lens vec
    vector<event_ind_t> v3;  // event classes
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

    GeneSegmentAlphabet gsa(VARIABLE, "testseg", alvec1, seqvec1);

    YMIR_ASSERT(gsa.name() == "testseg")
    YMIR_ASSERT(gsa.size() == 4)
    YMIR_ASSERT(gsa[0].sequence.size() == gsa["other"].sequence.size())
    YMIR_ASSERT(gsa[1].sequence == "ACT")
    YMIR_ASSERT(gsa[3].sequence == gsa["Vseg3"].sequence)

    gsa.appendPalindromicNucleotides(0, 0);
    YMIR_ASSERT2(gsa[1].sequence, "ACT")

    gsa.appendPalindromicNucleotides(2, 0);
    YMIR_ASSERT2(gsa[1].orig_sequence, "ACT")
    YMIR_ASSERT2(gsa[1].sequence, "GTACT")

    gsa.appendPalindromicNucleotides(0, 2);
    YMIR_ASSERT2(gsa[1].sequence, "ACTAG")

    gsa.appendPalindromicNucleotides(2, 2);
    YMIR_ASSERT2(gsa[2].orig_sequence, "GGG")
    YMIR_ASSERT2(gsa[2].sequence, "CCGGGCC")
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
    GeneSegmentAlphabet gsa_n(VARIABLE, "testseg", TEST_DATA_FOLDER + "RANDOM_FILE.txt", &ok);
    YMIR_ASSERT(!ok)
    GeneSegmentAlphabet gsa(VARIABLE, "testseg", TEST_DATA_FOLDER + "vgene.txt", &ok);
    YMIR_ASSERT(ok)

    YMIR_ASSERT(gsa.name() == "testseg")
    YMIR_ASSERT(gsa.size() == 4)
    YMIR_ASSERT(gsa[0].sequence.size() == gsa["other"].sequence.size())
    YMIR_ASSERT(gsa[1].sequence == "ACT")
    YMIR_ASSERT(gsa[3].sequence == gsa["Vseg3"].sequence)


    // assert write and than read again
    YMIR_ASSERT(gsa.write(TEST_DATA_FOLDER + "vgene_towrite.txt"))

    GeneSegmentAlphabet gsa1(VARIABLE, "testseg", TEST_DATA_FOLDER + "vgene_towrite.txt");

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
    seg_index_t *segs = new seg_index_t[10];
    segs[0] = 3;  // 3 Vs
    segs[1] = 2;  // 2 Js
    segs[2] = 2;  // 2 Ds

    segs[3] = 1;
    segs[4] = 2;
    segs[5] = 3;

    segs[6] = 4;
    segs[7] = 5;

    segs[8] = 6;
    segs[9] = 7;

    seq_len_t *alignments = new seq_len_t[27];
    // V
    alignments[0] = 8;
    alignments[1] = 9;
    alignments[2] = 5;

    alignments[3] = 11;
    alignments[4] = 12;
    alignments[5] = 3;

    alignments[6] = 14;
    alignments[7] = 15;
    alignments[8] = 6;

    // J
    alignments[9] = 17;
    alignments[10] = 18;
    alignments[11] = 7;

    alignments[12] = 20;
    alignments[13] = 21;
    alignments[14] = 3;

    // 2 X 3 tuples for D1
    alignments[15] = 31;
    alignments[16] = 32;
    alignments[17] = 4;

    alignments[18] = 41;
    alignments[19] = 42;
    alignments[20] = 7;

    // 2 X 3 tuples for D2
    alignments[21] = 51;
    alignments[22] = 52;
    alignments[23] = 6;

    alignments[24] = 60;
    alignments[25] = 61;
    alignments[26] = 4;

    seq_len_t *nd = new seq_len_t[2];
    nd[0] = 2;
    nd[1] = 2;

    Clonotype c("cloneseq", NUCLEOTIDE, VDJ_RECOMB, segs, alignments, nd);
    YMIR_ASSERT2(c.recombination(), VDJ_RECOMB)
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
    YMIR_ASSERT(c.numDivAlignments(0) == 2)

    YMIR_ASSERT2(c.getVarAlignment(0).gene_start(), 8)
    YMIR_ASSERT2(c.getVarAlignment(0).seq_end(), 13)
    YMIR_ASSERT2(c.getVarAlignment(0).length(), 5)

    YMIR_ASSERT2(c.getVarAlignment(2).gene_end(), 19)
    YMIR_ASSERT2(c.getVarAlignment(2).seq_start(), 15)
    YMIR_ASSERT2(c.getVarAlignment(2).length(), 6)

    YMIR_ASSERT2(c.getJoiAlignment(0).gene_start(), 17)
    YMIR_ASSERT2(c.getJoiAlignment(0).seq_end(), 24)
    YMIR_ASSERT2(c.getJoiAlignment(0).length(), 7)

    YMIR_ASSERT2(c.getJoiAlignment(1).gene_end(), 22)
    YMIR_ASSERT2(c.getJoiAlignment(1).seq_start(), 21)
    YMIR_ASSERT2(c.getJoiAlignment(1).length(), 3)
    
    YMIR_ASSERT2(c.getDivAlignment(0, 0).gene_start(), 31)
    YMIR_ASSERT2(c.getDivAlignment(0, 0).gene_end(), 34)
    YMIR_ASSERT2(c.getDivAlignment(0, 0).seq_start(), 32)
    YMIR_ASSERT2(c.getDivAlignment(0, 0).seq_end(), 35)

    YMIR_ASSERT2(c.getDivAlignment(0, 1).seq_start(), 42)
    YMIR_ASSERT2(c.getDivAlignment(0, 1).seq_end(), 48)

    YMIR_ASSERT2(c.getDivAlignment(1, 0).gene_start(), 51)
    YMIR_ASSERT2(c.getDivAlignment(1, 0).gene_end(), 56)
    YMIR_ASSERT2(c.getDivAlignment(1, 0).seq_start(), 52)
    YMIR_ASSERT2(c.getDivAlignment(1, 0).seq_end(), 57)

    YMIR_ASSERT2(c.getDivAlignment(1, 1).gene_start(), 60)
    YMIR_ASSERT2(c.getDivAlignment(1, 1).gene_end(), 63)
    YMIR_ASSERT2(c.getDivAlignment(1, 1).seq_start(), 61)
    YMIR_ASSERT2(c.getDivAlignment(1, 1).seq_end(), 64)
YMIR_TEST_END


YMIR_TEST_START(test_clonebuilder_clonealign)

    ClonotypeBuilder cb;

    cb.setNucleotideSeq();
    cb.setSequence("nuclseq");
    cb.addVarAlignment(10, 1, 3, 15)
            .addVarAlignment(11, 4, 1, 25)
            .addVarAlignment(12, 2, 6, 35)
            .addJoiAlignment(20, 1, 21, 10)
            .addDivAlignment(31, 8, 11, 2)
            .addDivAlignment(31, 8, 13, 2)
            .addDivAlignment(30, 1, 3, 2)
            .addDivAlignment(30, 1, 5, 2)
            .addDivAlignment(32, 1, 5, 2);

    cb.setRecombination(VDJ_RECOMB);
    Clonotype c = cb.buildClonotype();

    YMIR_ASSERT2(c.sequence(), "nuclseq")
    YMIR_ASSERT2(c.sequence_type(), NUCLEOTIDE)
    YMIR_ASSERT2(c.recombination(), VDJ_RECOMB)

    YMIR_ASSERT2(c.getVar(0), 10)
    YMIR_ASSERT(c.getVarAlignment(0) == Alignment(1, 3, 15))

    YMIR_ASSERT2(c.getVar(1), 11)
    YMIR_ASSERT(c.getVarAlignment(1) == Alignment(4, 1, 25))

    YMIR_ASSERT2(c.getVar(2), 12)
    YMIR_ASSERT(c.getVarAlignment(2) == Alignment(2, 6, 35))

    YMIR_ASSERT2(c.getJoi(0), 20)
    YMIR_ASSERT(c.getJoiAlignment(0) == Alignment(1, 21, 10))

    YMIR_ASSERT2(c.getDiv(0), 31)
    YMIR_ASSERT(c.getDivAlignment(0, 0) == Alignment(8, 11, 2))
    YMIR_ASSERT(c.getDivAlignment(0, 1) == Alignment(8, 13, 2))
    YMIR_ASSERT2(c.getDiv(1), 30)
    YMIR_ASSERT(c.getDivAlignment(1, 0) == Alignment(1, 3, 2))
    YMIR_ASSERT(c.getDivAlignment(1, 1) == Alignment(1, 5, 2))
    YMIR_ASSERT2(c.getDiv(2), 32)
    YMIR_ASSERT(c.getDivAlignment(2, 0) == Alignment(1, 5, 2))
    YMIR_ASSERT2(c.nVar(), 3)
    YMIR_ASSERT2(c.nJoi(), 1)
    YMIR_ASSERT2(c.nDiv(), 3)

YMIR_TEST_END


YMIR_TEST_START(test_nogap_alignment_vector_no_err)

    NoGapAlignmentVector vec;

    YMIR_ASSERT2(vec.size(), 0)

    vec.addAlignment(11, 1, 2, 3);

    YMIR_ASSERT2(vec.size(), 1)
    YMIR_ASSERT2(vec.pattern_start(0), 1)
    YMIR_ASSERT2(vec.text_start(0), 2)
    YMIR_ASSERT2(vec.len(0), 3)
    YMIR_ASSERT2(vec.id(0), 11)

    vec.addAlignment(14, 4, 5, 6);

    YMIR_ASSERT2(vec.size(), 2)
    YMIR_ASSERT2(vec.pattern_start(0), 1)
    YMIR_ASSERT2(vec.text_start(0), 2)
    YMIR_ASSERT2(vec.len(0), 3)
    YMIR_ASSERT2(vec.id(0), 11)

    YMIR_ASSERT2(vec.pattern_start(1), 4)
    YMIR_ASSERT2(vec.text_start(1), 5)
    YMIR_ASSERT2(vec.len(1), 6)
    YMIR_ASSERT2(vec.id(1), 14)

    vec.addAlignment(18, 8, 9, 10);

    YMIR_ASSERT2(vec.size(), 3)
    YMIR_ASSERT2(vec.pattern_start(1), 4)
    YMIR_ASSERT2(vec.text_start(1), 5)
    YMIR_ASSERT2(vec.len(1), 6)
    YMIR_ASSERT2(vec.id(1), 14)

    YMIR_ASSERT2(vec.pattern_start(2), 8)
    YMIR_ASSERT2(vec.text_start(2), 9)
    YMIR_ASSERT2(vec.len(2), 10)
    YMIR_ASSERT2(vec.id(2), 18)

YMIR_TEST_END


YMIR_TEST_START(test_nogap_alignment_vector_errors)

    NoGapAlignmentVector vec;

    YMIR_ASSERT2(vec.size(), 0)

    AlignmentVectorBase::events_storage_t events1 {false, true, true};

    vec.addAlignment(11, 1, 2, events1);
    YMIR_ASSERT2(vec.size(), 1)
    YMIR_ASSERT2(vec.pattern_start(0), 1)
    YMIR_ASSERT2(vec.text_start(0), 2)
    YMIR_ASSERT2(vec.len(0), 3)
    YMIR_ASSERT2(vec.id(0), 11)
    YMIR_ASSERT2(vec.isMismatch(0, 1), false)
    YMIR_ASSERT2(vec.isMismatch(0, 2), true)
    YMIR_ASSERT2(vec.isMismatch(0, 3), true)

    AlignmentVectorBase::events_storage_t events2 {false, false, true, false};

    vec.addAlignment(13, 3, 4, events2);
    YMIR_ASSERT2(vec.size(), 2)
    YMIR_ASSERT2(vec.pattern_start(0), 1)
    YMIR_ASSERT2(vec.text_start(0), 2)
    YMIR_ASSERT2(vec.len(0), 3)
    YMIR_ASSERT2(vec.id(0), 11)
    YMIR_ASSERT2(vec.isMismatch(0, 1), false)
    YMIR_ASSERT2(vec.isMismatch(0, 2), true)
    YMIR_ASSERT2(vec.isMismatch(0, 3), true)

    YMIR_ASSERT2(vec.pattern_start(1), 3)
    YMIR_ASSERT2(vec.text_start(1), 4)
    YMIR_ASSERT2(vec.len(1), 4)
    YMIR_ASSERT2(vec.id(1), 13)
    YMIR_ASSERT2(vec.isMismatch(1, 1), false)
    YMIR_ASSERT2(vec.isMismatch(1, 2), false)
    YMIR_ASSERT2(vec.isMismatch(1, 3), true)
    YMIR_ASSERT2(vec.isMismatch(1, 4), false)


    NoGapAlignmentVector vec2;

    AlignmentVectorBase::events_storage_t events3 {true, true, true, false, true, false};
    vec2.addAlignment(21, 5, 6, events3);

    vec.extend(vec2);
    YMIR_ASSERT2(vec.size(), 3)

    YMIR_ASSERT2(vec.pattern_start(0), 1)
    YMIR_ASSERT2(vec.text_start(0), 2)
    YMIR_ASSERT2(vec.len(0), 3)
    YMIR_ASSERT2(vec.id(0), 11)
    YMIR_ASSERT2(vec.isMismatch(0, 1), false)
    YMIR_ASSERT2(vec.isMismatch(0, 2), true)
    YMIR_ASSERT2(vec.isMismatch(0, 3), true)

    YMIR_ASSERT2(vec.pattern_start(1), 3)
    YMIR_ASSERT2(vec.text_start(1), 4)
    YMIR_ASSERT2(vec.len(1), 4)
    YMIR_ASSERT2(vec.id(1), 13)
    YMIR_ASSERT2(vec.isMismatch(1, 1), false)
    YMIR_ASSERT2(vec.isMismatch(1, 2), false)
    YMIR_ASSERT2(vec.isMismatch(1, 3), true)
    YMIR_ASSERT2(vec.isMismatch(1, 4), false)

    YMIR_ASSERT2(vec.pattern_start(2), 5)
    YMIR_ASSERT2(vec.text_start(2), 6)
    YMIR_ASSERT2(vec.len(2), 6)
    YMIR_ASSERT2(vec.id(2), 21)
    YMIR_ASSERT2(vec.isMismatch(2, 1), true)
    YMIR_ASSERT2(vec.isMismatch(2, 2), true)
    YMIR_ASSERT2(vec.isMismatch(2, 3), true)
    YMIR_ASSERT2(vec.isMismatch(2, 4), false)
    YMIR_ASSERT2(vec.isMismatch(2, 5), true)
    YMIR_ASSERT2(vec.isMismatch(2, 6), false)

YMIR_TEST_END


YMIR_TEST_START(test_gapped_alignment_vector)

    GappedAlignmentVector vec;

    YMIR_ASSERT2(vec.size(), 0)

    // match mismatch ins ins
    AlignmentVectorBase::events_storage_t events1;
    add_match(&events1);
    add_mismatch(&events1);
    add_ins(&events1);
    add_ins(&events1);

    vec.addAlignment(11, 1, 2, events1);
    YMIR_ASSERT2(vec.size(), 1)
    YMIR_ASSERT2(vec.pattern_start(0), 1)
    YMIR_ASSERT2(vec.text_start(0), 2)
    YMIR_ASSERT2(vec.len(0), 4)
    YMIR_ASSERT2(vec.id(0), 11)
    YMIR_ASSERT(vec.isMatch(0, 1))
    YMIR_ASSERT(vec.isMismatch(0, 2))
    YMIR_ASSERT(vec.isIns(0, 3))
    YMIR_ASSERT(vec.isIns(0, 4))

    // mismath ins match del del
    AlignmentVectorBase::events_storage_t events2;
    add_mismatch(&events2);
    add_ins(&events2);
    add_match(&events2);
    add_del(&events2);
    add_del(&events2);

    vec.addAlignment(14, 3, 4, events2);
    YMIR_ASSERT2(vec.size(), 2)
    YMIR_ASSERT2(vec.pattern_start(0), 1)
    YMIR_ASSERT2(vec.text_start(0), 2)
    YMIR_ASSERT2(vec.len(0), 4)
    YMIR_ASSERT2(vec.id(0), 11)
    YMIR_ASSERT(vec.isMatch(0, 1))
    YMIR_ASSERT(vec.isMismatch(0, 2))
    YMIR_ASSERT(vec.isIns(0, 3))
    YMIR_ASSERT(vec.isIns(0, 4))

    YMIR_ASSERT2(vec.pattern_start(1), 3)
    YMIR_ASSERT2(vec.text_start(1), 4)
    YMIR_ASSERT2(vec.len(1), 5)
    YMIR_ASSERT2(vec.id(1), 14)
    YMIR_ASSERT(vec.isMismatch(1, 1))
    YMIR_ASSERT(vec.isIns(1, 2))
    YMIR_ASSERT(vec.isMatch(1, 3))
    YMIR_ASSERT(vec.isDel(1, 4))
    YMIR_ASSERT(vec.isDel(1, 5))


    GappedAlignmentVector vec2;
    AlignmentVectorBase::events_storage_t events3;
    add_ins(&events3);
    add_del(&events3);
    add_match(&events3);
    add_mismatch(&events3);
    add_del(&events3);
    add_ins(&events3);
    vec2.addAlignment(21, 5, 6, events3);

    vec.extend(vec2);

    YMIR_ASSERT2(vec.size(), 3)

    YMIR_ASSERT2(vec.pattern_start(0), 1)
    YMIR_ASSERT2(vec.text_start(0), 2)
    YMIR_ASSERT2(vec.len(0), 4)
    YMIR_ASSERT2(vec.id(0), 11)
    YMIR_ASSERT(vec.isMatch(0, 1))
    YMIR_ASSERT(vec.isMismatch(0, 2))
    YMIR_ASSERT(vec.isIns(0, 3))
    YMIR_ASSERT(vec.isIns(0, 4))

    YMIR_ASSERT2(vec.pattern_start(1), 3)
    YMIR_ASSERT2(vec.text_start(1), 4)
    YMIR_ASSERT2(vec.len(1), 5)
    YMIR_ASSERT2(vec.id(1), 14)
    YMIR_ASSERT(vec.isMismatch(1, 1))
    YMIR_ASSERT(vec.isIns(1, 2))
    YMIR_ASSERT(vec.isMatch(1, 3))
    YMIR_ASSERT(vec.isDel(1, 4))
    YMIR_ASSERT(vec.isDel(1, 5))

    YMIR_ASSERT2(vec.pattern_start(2), 5)
    YMIR_ASSERT2(vec.text_start(2), 6)
    YMIR_ASSERT2(vec.len(2), 6)
    YMIR_ASSERT2(vec.id(2), 21)    
    YMIR_ASSERT(vec.isIns(2, 1))
    YMIR_ASSERT(vec.isDel(2, 2))
    YMIR_ASSERT(vec.isMatch(2, 3))
    YMIR_ASSERT(vec.isMismatch(2, 4))
    YMIR_ASSERT(vec.isDel(2, 5))
    YMIR_ASSERT(vec.isIns(2, 6))

YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_simple_vj)

    VDJAlignmentBuilder builder;

    builder.addVarAlignment(11, 1, 2, 3)
           .addVarAlignment(12, 4, 5, 6)
           .addJoiAlignment(31, 7, 8, 9)
           .addJoiAlignment(33, 10, 11, 12)
           .addJoiAlignment(35, 13, 14, 15);

    VDJAlignment algn = builder.build();

    YMIR_ASSERT2(algn.nVar(), 2)
    YMIR_ASSERT2(algn.nJoi(), 3)
    YMIR_ASSERT2(algn.nDiv(), 0)

    YMIR_ASSERT2(algn.getVar(0), 11)
    YMIR_ASSERT2(algn.getVar(1), 12)
    YMIR_ASSERT2(algn.getJoi(0), 31)  
    YMIR_ASSERT2(algn.getJoi(1), 33)  
    YMIR_ASSERT2(algn.getJoi(2), 35)

    YMIR_ASSERT2(algn.getVarGeneStart(0), 1)
    YMIR_ASSERT2(algn.getVarSeqStart(0), 2)
    YMIR_ASSERT2(algn.getVarLen(0), 3)
    YMIR_ASSERT2(algn.getVarGeneStart(1), 4)
    YMIR_ASSERT2(algn.getVarSeqStart(1), 5)
    YMIR_ASSERT2(algn.getVarLen(1), 6)

    YMIR_ASSERT2(algn.getJoiGeneStart(0), 7)
    YMIR_ASSERT2(algn.getJoiSeqStart(0), 8)
    YMIR_ASSERT2(algn.getJoiLen(0), 9)
    YMIR_ASSERT2(algn.getJoiGeneStart(2), 13)
    YMIR_ASSERT2(algn.getJoiSeqStart(2), 14)
    YMIR_ASSERT2(algn.getJoiLen(2), 15)


    builder.addVarAlignment(41, 41, 42, 43)
           .addVarAlignment(42, 44, 45, 46)
           .addVarAlignment(49, 57, 58, 59)
           .addJoiAlignment(53, 1, 2, 3)
           .addJoiAlignment(55, 4, 5, 6);

    algn = builder.build();

    YMIR_ASSERT2(algn.nVar(), 3)
    YMIR_ASSERT2(algn.nJoi(), 2)
    YMIR_ASSERT2(algn.nDiv(), 0)

    YMIR_ASSERT2(algn.getVar(0), 41)
    YMIR_ASSERT2(algn.getVar(1), 42)
    YMIR_ASSERT2(algn.getVar(2), 49)  
    YMIR_ASSERT2(algn.getJoi(0), 53)  
    YMIR_ASSERT2(algn.getJoi(1), 55)

    YMIR_ASSERT2(algn.getVarGeneStart(0), 41)
    YMIR_ASSERT2(algn.getVarSeqStart(0), 42)
    YMIR_ASSERT2(algn.getVarLen(0), 43)
    YMIR_ASSERT2(algn.getVarGeneStart(1), 44)
    YMIR_ASSERT2(algn.getVarSeqStart(1), 45)
    YMIR_ASSERT2(algn.getVarLen(1), 46)
    YMIR_ASSERT2(algn.getVarGeneStart(2), 57)
    YMIR_ASSERT2(algn.getVarSeqStart(2), 58)
    YMIR_ASSERT2(algn.getVarLen(2), 59)

    YMIR_ASSERT2(algn.getJoiGeneStart(0), 1)
    YMIR_ASSERT2(algn.getJoiSeqStart(0), 2)
    YMIR_ASSERT2(algn.getJoiLen(0), 3)
    YMIR_ASSERT2(algn.getJoiGeneStart(1), 4)
    YMIR_ASSERT2(algn.getJoiSeqStart(1), 5)
    YMIR_ASSERT2(algn.getJoiLen(1), 6)

YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_simple_vdj)

    VDJAlignmentBuilder builder;

    // VDJ

    builder.addVarAlignment(11, 1, 2, 3)
           .addVarAlignment(12, 4, 5, 6)
           .addJoiAlignment(31, 7, 8, 9)
           .addJoiAlignment(33, 10, 11, 12)
           .addJoiAlignment(35, 13, 14, 15)
           .addDivAlignment(41, 20, 21, 22)
           .addDivAlignment(41, 23, 24, 25)
           .addDivAlignment(45, 30, 31, 32)
           .addDivAlignment(45, 33, 34, 35)
           .addDivAlignment(45, 36, 37, 38);

    VDJAlignment algn = builder.build();

    YMIR_ASSERT2(algn.nVar(), 2)
    YMIR_ASSERT2(algn.nJoi(), 3)
    YMIR_ASSERT2(algn.nDiv(), 2)
    YMIR_ASSERT2(algn.numDivAlignments(0), 2)
    YMIR_ASSERT2(algn.numDivAlignments(1), 3)

    YMIR_ASSERT2(algn.getVar(0), 11)
    YMIR_ASSERT2(algn.getVar(1), 12)
    YMIR_ASSERT2(algn.getJoi(0), 31)
    YMIR_ASSERT2(algn.getJoi(1), 33)
    YMIR_ASSERT2(algn.getJoi(2), 35)
    YMIR_ASSERT2(algn.getDiv(0), 41)
    YMIR_ASSERT2(algn.getDiv(1), 45)

    YMIR_ASSERT2(algn.getVarGeneStart(0), 1)
    YMIR_ASSERT2(algn.getVarSeqStart(0), 2)
    YMIR_ASSERT2(algn.getVarLen(0), 3)
    YMIR_ASSERT2(algn.getVarGeneStart(1), 4)
    YMIR_ASSERT2(algn.getVarSeqStart(1), 5)
    YMIR_ASSERT2(algn.getVarLen(1), 6)

    YMIR_ASSERT2(algn.getJoiGeneStart(0), 7)
    YMIR_ASSERT2(algn.getJoiSeqStart(0), 8)
    YMIR_ASSERT2(algn.getJoiLen(0), 9)
    YMIR_ASSERT2(algn.getJoiGeneStart(2), 13)
    YMIR_ASSERT2(algn.getJoiSeqStart(2), 14)
    YMIR_ASSERT2(algn.getJoiLen(2), 15)

    YMIR_ASSERT2(algn.getDivGeneStart(0, 0), 20)
    YMIR_ASSERT2(algn.getDivSeqStart(0, 0), 21)
    YMIR_ASSERT2(algn.getDivLen(0, 0), 22)
    YMIR_ASSERT2(algn.getDivGeneStart(0, 1), 23)
    YMIR_ASSERT2(algn.getDivSeqStart(0, 1), 24)
    YMIR_ASSERT2(algn.getDivLen(0, 1), 25)

    YMIR_ASSERT2(algn.getDivGeneStart(1, 0), 30)
    YMIR_ASSERT2(algn.getDivSeqStart(1, 0), 31)
    YMIR_ASSERT2(algn.getDivLen(1, 0), 32)
    YMIR_ASSERT2(algn.getDivGeneStart(1, 2), 36)
    YMIR_ASSERT2(algn.getDivSeqStart(1, 2), 37)
    YMIR_ASSERT2(algn.getDivLen(1, 2), 38)

YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_vector_vj)

    VDJAlignmentBuilder builder;

    // V

    NoGapAlignmentVector vec1;
    AlignmentVectorBase::events_storage_t events11 {false, true, true};
    vec1.addAlignment(11, 1, 2, events11);
    AlignmentVectorBase::events_storage_t events12 {false, false, true, 
                                                    true, true, false};
    vec1.addAlignment(12, 4, 5, events12);


    // J

    NoGapAlignmentVector vec2;
    AlignmentVectorBase::events_storage_t events21 {false, true, true, 
                                                    false, true, true, 
                                                    false, true, true};
    vec2.addAlignment(31, 7, 8, events21);
    AlignmentVectorBase::events_storage_t events22 {true, false, true, 
                                                    true, true, false, 
                                                    false, false, true, 
                                                    true, true, false};
    vec2.addAlignment(33, 10, 11, events22);

    NoGapAlignmentVector vec3;
    AlignmentVectorBase::events_storage_t events31 {false, true, true, 
                                                    false, true, true, 
                                                    false, true, true,
                                                    false, true, true,
                                                    false, true, true};
    vec3.addAlignment(35, 13, 14, events31);


    builder.addVarAlignment(vec1)
           .addJoiAlignment(vec2)
           .addJoiAlignment(vec3);

    VDJAlignment algn = builder.build();

    YMIR_ASSERT2(algn.nVar(), 2)
    YMIR_ASSERT2(algn.nJoi(), 3)
    YMIR_ASSERT2(algn.nDiv(), 0)

    YMIR_ASSERT2(algn.getVar(0), 11)
    YMIR_ASSERT2(algn.getVar(1), 12)
    YMIR_ASSERT2(algn.getJoi(0), 31)
    YMIR_ASSERT2(algn.getJoi(1), 33)
    YMIR_ASSERT2(algn.getJoi(2), 35)

    YMIR_ASSERT2(algn.getVarGeneStart(0), 1)
    YMIR_ASSERT2(algn.getVarSeqStart(0), 2)
    YMIR_ASSERT2(algn.getVarLen(0), 3)
    YMIR_ASSERT(!algn.isVarMismatch(0, 1))
    YMIR_ASSERT(algn.isVarMismatch(0, 2))
    YMIR_ASSERT(algn.isVarMismatch(0, 3))

    YMIR_ASSERT2(algn.getVarGeneStart(1), 4)
    YMIR_ASSERT2(algn.getVarSeqStart(1), 5)
    YMIR_ASSERT2(algn.getVarLen(1), 6)
    YMIR_ASSERT(!algn.isVarMismatch(1, 1))
    YMIR_ASSERT(algn.isVarMismatch(1, 3))
    YMIR_ASSERT(!algn.isVarMismatch(1, 6))

    YMIR_ASSERT2(algn.getJoiGeneStart(0), 7)
    YMIR_ASSERT2(algn.getJoiSeqStart(0), 8)
    YMIR_ASSERT2(algn.getJoiLen(0), 9)
    YMIR_ASSERT(!algn.isJoiMismatch(0, 1))
    YMIR_ASSERT(algn.isJoiMismatch(0, 5))
    YMIR_ASSERT(algn.isJoiMismatch(0, 9))

    YMIR_ASSERT(algn.isJoiMismatch(1, 1))
    YMIR_ASSERT(algn.isJoiMismatch(1, 4))
    YMIR_ASSERT(!algn.isJoiMismatch(1, 12))

    YMIR_ASSERT2(algn.getJoiGeneStart(2), 13)
    YMIR_ASSERT2(algn.getJoiSeqStart(2), 14)
    YMIR_ASSERT2(algn.getJoiLen(2), 15)
    YMIR_ASSERT(!algn.isJoiMismatch(2, 1))
    YMIR_ASSERT(!algn.isJoiMismatch(2, 4))

    YMIR_ASSERT(algn.isJoiMismatch(2, 14))
    YMIR_ASSERT(algn.isJoiMismatch(2, 15))

YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_vector_vdj)

    VDJAlignmentBuilder builder;

    // V

    NoGapAlignmentVector vec1;
    AlignmentVectorBase::events_storage_t events11 {false, true, true};
    vec1.addAlignment(11, 1, 2, events11);
    AlignmentVectorBase::events_storage_t events12 {false, false, true, 
                                                    true, true, false};
    vec1.addAlignment(12, 4, 5, events12);


    // J

    NoGapAlignmentVector vec2;
    AlignmentVectorBase::events_storage_t events21 {false, true, true, 
                                                    false, true, true, 
                                                    false, true, true};
    vec2.addAlignment(31, 7, 8, events21);
    AlignmentVectorBase::events_storage_t events22 {true, false, true, 
                                                    true, true, false, 
                                                    false, false, true, 
                                                    true, true, false};
    vec2.addAlignment(33, 10, 11, events22);

    NoGapAlignmentVector vec3;
    AlignmentVectorBase::events_storage_t events31 {false, true, true, 
                                                    false, true, true, 
                                                    false, true, true,
                                                    false, true, true,
                                                    false, true, true};
    vec3.addAlignment(35, 13, 14, events31);


    // D

    NoGapAlignmentVector vec4;
    AlignmentVectorBase::events_storage_t events41 {false, true, true, 
                                                    false, true, true};
    AlignmentVectorBase::events_storage_t events42 {false, true, false};
    vec4.addAlignment(43, 13, 14, events41);
    vec4.addAlignment(43, 15, 16, events42);

    NoGapAlignmentVector vec5;
    AlignmentVectorBase::events_storage_t events51 {false, false, true, true, 
                                                    false, true, true};
    vec5.addAlignment(45, 17, 18, events51);


    builder.addVarAlignment(vec1)
           .addJoiAlignment(vec2)
           .addJoiAlignment(vec3)
           .addDivAlignment(vec4)
           .addDivAlignment(vec5);

    VDJAlignment algn = builder.build();

    YMIR_ASSERT2(algn.nVar(), 2)
    YMIR_ASSERT2(algn.nJoi(), 3)
    YMIR_ASSERT2(algn.nDiv(), 2)

    YMIR_ASSERT2(algn.numDivAlignments(0), 2)
    YMIR_ASSERT2(algn.numDivAlignments(1), 1)

    YMIR_ASSERT2(algn.getVar(0), 11)
    YMIR_ASSERT2(algn.getVar(1), 12)
    YMIR_ASSERT2(algn.getJoi(0), 31)
    YMIR_ASSERT2(algn.getJoi(1), 33)
    YMIR_ASSERT2(algn.getJoi(2), 35)

    YMIR_ASSERT2(algn.getVarGeneStart(0), 1)
    YMIR_ASSERT2(algn.getVarSeqStart(0), 2)
    YMIR_ASSERT2(algn.getVarLen(0), 3)
    YMIR_ASSERT(!algn.isVarMismatch(0, 1))
    YMIR_ASSERT(algn.isVarMismatch(0, 2))
    YMIR_ASSERT(algn.isVarMismatch(0, 3))

    YMIR_ASSERT2(algn.getVarGeneStart(1), 4)
    YMIR_ASSERT2(algn.getVarSeqStart(1), 5)
    YMIR_ASSERT2(algn.getVarLen(1), 6)
    YMIR_ASSERT(!algn.isVarMismatch(1, 1))
    YMIR_ASSERT(algn.isVarMismatch(1, 3))
    YMIR_ASSERT(!algn.isVarMismatch(1, 6))

    YMIR_ASSERT2(algn.getJoiGeneStart(0), 7)
    YMIR_ASSERT2(algn.getJoiSeqStart(0), 8)
    YMIR_ASSERT2(algn.getJoiLen(0), 9)
    YMIR_ASSERT(!algn.isJoiMismatch(0, 1))
    YMIR_ASSERT(algn.isJoiMismatch(0, 5))
    YMIR_ASSERT(algn.isJoiMismatch(0, 9))

    YMIR_ASSERT(algn.isJoiMismatch(1, 1))
    YMIR_ASSERT(algn.isJoiMismatch(1, 4))
    YMIR_ASSERT(!algn.isJoiMismatch(1, 12))

    YMIR_ASSERT2(algn.getJoiGeneStart(2), 13)
    YMIR_ASSERT2(algn.getJoiSeqStart(2), 14)
    YMIR_ASSERT2(algn.getJoiLen(2), 15)
    YMIR_ASSERT(!algn.isJoiMismatch(2, 1))
    YMIR_ASSERT(!algn.isJoiMismatch(2, 4))

    YMIR_ASSERT(algn.isJoiMismatch(2, 14))
    YMIR_ASSERT(algn.isJoiMismatch(2, 15))

    // YMIR_ASSERT2(algn.getDivGeneStart(0, 0), 13)
    // YMIR_ASSERT2(algn.getDivSeqStart(0, 0), 14)
    // YMIR_ASSERT2(algn.getDivLen(0, 0), 6)
    // YMIR_ASSERT(!algn.isDivMismatch(0, 0, 1))
    // YMIR_ASSERT(algn.isDivMismatch(0, 0, 2))
    // YMIR_ASSERT(algn.isDivMismatch(0, 0, 3))
    // YMIR_ASSERT(!algn.isDivMismatch(0, 0, 4))
    // YMIR_ASSERT(algn.isDivMismatch(0, 0, 5))
    // YMIR_ASSERT(algn.isDivMismatch(0, 0, 6))

    // YMIR_ASSERT2(algn.getDivGeneStart(0, 1), 15)
    // YMIR_ASSERT2(algn.getDivSeqStart(0, 1), 16)
    // YMIR_ASSERT2(algn.getDivLen(0, 1), 3)
    // YMIR_ASSERT(!algn.isDivMismatch(0, 1, 1))
    // YMIR_ASSERT(algn.isDivMismatch(0, 1, 2))
    // YMIR_ASSERT(!algn.isDivMismatch(0, 1, 3))

    // YMIR_ASSERT2(algn.getDivGeneStart(1, 0), 17)
    // YMIR_ASSERT2(algn.getDivSeqStart(1, 0), 18)
    // YMIR_ASSERT2(algn.getDivLen(1, 0), 7)
    // YMIR_ASSERT(!algn.isDivMismatch(1, 0, 1))
    // YMIR_ASSERT(!algn.isDivMismatch(1, 0, 2))
    // YMIR_ASSERT(algn.isDivMismatch(1, 0, 3))
    // YMIR_ASSERT(algn.isDivMismatch(1, 0, 4))
    // YMIR_ASSERT(!algn.isDivMismatch(1, 0, 5))
    // YMIR_ASSERT(algn.isDivMismatch(1, 0, 6))
    // YMIR_ASSERT(algn.isDivMismatch(1, 0, 7))

YMIR_TEST_END


YMIR_TEST_START(test_naive_cdr3_nuc_aligner)

    NoGapAlignmentVector vec;

    vector<string> avec1 {"V1", "V2", "V3", "V4"};
    vector<string> svec1 {"ACGTT", "ACGT", "ACG", "TTT"};

    vector<string> avec2 {"J1", "J2", "J3", "J4"};
    vector<string> svec2 {"CGT", "TACGT", "TTCGT", "TTTTT"};

    vector<string> avec3 {"D1", "D2", "D3"};
    vector<string> svec3 {"AA", "AACCTT", "ACT"};

    VDJRecombinationGenes genes("V", avec1, svec1, "J", avec2, svec2, "D", avec3, svec3);

    NaiveCDR3NucleotideAligner nna(genes, VDJAlignerParameters(1, 3));

    YMIR_ASSERT2(nna.alignVar(1, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(1, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(1, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignVar(2, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(2, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(2, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignVar(3, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(3, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(3, "ACGT").len(0), 3)

    YMIR_ASSERT2(nna.alignVar(4, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(4, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(4, "ACGT").len(0), 0)

    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").pattern_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").len(0), 3)

    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").text_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").pattern_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").text_start(0), 3)
    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").len(0), 3)

    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").pattern_start(0), 4)
    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").text_start(0), 5)
    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").len(0), 1)

    YMIR_ASSERT2(nna.alignJoi(4, "ACGG").len(0), 0)

    YMIR_ASSERT2(nna.alignDiv(1, "TTAATAA").size(), 0)

    NaiveCDR3NucleotideAligner nna2(genes, VDJAlignerParameters(1, 2));
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").size(), 2)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").text_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").pattern_start(0), 3)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").len(0), 2)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").text_start(1), 1)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").pattern_start(1), 6)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").len(1), 2)

    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").size(), 3)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(0), 2)

    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(1), 5)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(1), 5)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(1), 2)

    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(2), 5)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(2), 12)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(2), 2)

    YMIR_ASSERT2(nna2.alignDiv(3, "ACTGACGACGGTATCTAC").size(), 5)

YMIR_TEST_END


YMIR_TEST_START(test_cdr3_nuc_aligner)

    NoGapAlignmentVector vec;

    vector<string> avec1 {"V1", "V2", "V3", "V4"};
    vector<string> svec1 {"ACGTT", "ACGT", "ACG", "TTT"};

    vector<string> avec2 {"J1", "J2", "J3", "J4"};
    vector<string> svec2 {"CGT", "TACGT", "TTCGT", "TTTTT"};

    vector<string> avec3 {"D1", "D2", "D3"};
    vector<string> svec3 {"AA", "AACCTT", "ACT"};

    VDJRecombinationGenes genes("V", avec1, svec1, "J", avec2, svec2, "D", avec3, svec3);

    CDR3NucleotideAligner nna(genes, VDJAlignerParameters(1, 3));

    YMIR_ASSERT2(nna.alignVar(1, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(1, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(1, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignVar(2, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(2, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(2, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignVar(3, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(3, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(3, "ACGT").len(0), 3)

    YMIR_ASSERT2(nna.alignVar(4, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(4, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(4, "ACGT").len(0), 3)
    YMIR_ASSERT(nna.alignVar(4, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(nna.alignVar(4, "ACGT").isMismatch(0, 2))
    YMIR_ASSERT(nna.alignVar(4, "ACGT").isMismatch(0, 3))

    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").pattern_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").len(0), 3)
    YMIR_ASSERT(!nna.alignJoi(1, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(!nna.alignJoi(1, "ACGT").isMismatch(0, 2))
    YMIR_ASSERT(!nna.alignJoi(1, "ACGT").isMismatch(0, 3))

    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").text_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").len(0), 4)
    YMIR_ASSERT(!nna.alignJoi(2, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(!nna.alignJoi(2, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(!nna.alignJoi(2, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(!nna.alignJoi(2, "ACGT").isMismatch(0, 1))

    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").text_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").len(0), 4)
    YMIR_ASSERT(nna.alignJoi(3, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(!nna.alignJoi(3, "ACGT").isMismatch(0, 2))
    YMIR_ASSERT(!nna.alignJoi(3, "ACGT").isMismatch(0, 3))
    YMIR_ASSERT(!nna.alignJoi(3, "ACGT").isMismatch(0, 4))

    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").text_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").len(0), 4)
    YMIR_ASSERT(nna.alignJoi(4, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(nna.alignJoi(4, "ACGT").isMismatch(0, 2))
    YMIR_ASSERT(nna.alignJoi(4, "ACGT").isMismatch(0, 3))
    YMIR_ASSERT(!nna.alignJoi(4, "ACGT").isMismatch(0, 4))

    YMIR_ASSERT2(nna.alignJoi(4, "ACGG").len(0), 4)
    YMIR_ASSERT(nna.alignJoi(4, "ACGG").isMismatch(0, 1))
    YMIR_ASSERT(nna.alignJoi(4, "ACGG").isMismatch(0, 2))
    YMIR_ASSERT(nna.alignJoi(4, "ACGG").isMismatch(0, 3))
    YMIR_ASSERT(nna.alignJoi(4, "ACGG").isMismatch(0, 4))

    // YMIR_ASSERT2(nna.alignDiv(1, "TTAATAA").size(), 0)

    // NaiveCDR3NucleotideAligner nna2(genes, VDJAlignerParameters(1, 2));
    // YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").size(), 2)
    // YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").text_start(0), 1)
    // YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").pattern_start(0), 3)
    // YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").len(0), 2)
    // YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").text_start(1), 1)
    // YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").pattern_start(1), 6)
    // YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").len(1), 2)

    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").size(), 3)
    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(0), 1)
    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(0), 1)
    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(0), 2)

    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(1), 5)
    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(1), 5)
    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(1), 2)

    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(2), 5)
    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(2), 12)
    // YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(2), 2)

    // YMIR_ASSERT2(nna2.alignDiv(3, "ACTGACGACGGTATCTAC").size(), 5)

YMIR_TEST_END


YMIR_TEST_START(test_cdr3_aa_aligner)

    YMIR_ASSERT(false)

    // NaiveAminoAcidAligner naa;

    // // {'S', "TCT"}, {'S', "TCC"}, {'S', "TCA"}, {'S', "TCG"}, {'S', "AGT"}, {'S', "AGC"},
    // // {'R', "CGT"}, {'R', "CGC"}, {'R', "CGA"}, {'R', "CGG"}, {'R', "AGA"}, {'R', "AGG"},
    // YMIR_ASSERT2(naa.align5end("GGG", "SR"), 0)
    // YMIR_ASSERT2(naa.align5end("TC", "SR"), 2)
    // YMIR_ASSERT2(naa.align5end("TCA", "SR"), 3)
    // YMIR_ASSERT2(naa.align5end("TCGCG", "SR"), 5)

    // YMIR_ASSERT2(naa.align5end("AG", "SR"), 2)
    // YMIR_ASSERT2(naa.align5end("AGTAGA", "SR"), 6)
    // YMIR_ASSERT2(naa.align5end("AGTAGG", "SR"), 6)

    // YMIR_ASSERT2(naa.align3end("T", "SR"), 1)
    // YMIR_ASSERT2(naa.align3end("CT", "SR"), 1)
    // YMIR_ASSERT2(naa.align3end("GT", "SR"), 2)

    // YMIR_ASSERT2(naa.align3end("AAA", "SR"), 1)
    // YMIR_ASSERT2(naa.align3end("TGA", "SR"), 2)
    // YMIR_ASSERT2(naa.align3end("TTTAGG", "SR"), 4)
    // YMIR_ASSERT2(naa.align3end("AGTAGG", "SR"), 6)

    // YMIR_ASSERT(naa.alignLocal("TGATGAA", "SR").size() != 0)

YMIR_TEST_END


YMIR_TEST_START(test_sw_aligner)

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_swng_aligner)

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_errcorr_aligner)

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_writer)

    // RepertoireWriter writer;

    // vector<string> alvec1;
    // vector<string> seqvec1;
    // alvec1.push_back("Vseg1");
    // alvec1.push_back("Vseg2");
    // alvec1.push_back("Vseg3");
    // seqvec1.push_back("CCCG");
    // seqvec1.push_back("GGG");
    // seqvec1.push_back("AGGCGAG");

    // vector<string> alvec2;
    // vector<string> seqvec2;
    // alvec2.push_back("Jseg1");
    // alvec2.push_back("Jseg2");
    // alvec2.push_back("Jseg3");
    // seqvec2.push_back("CCGTTT");
    // seqvec2.push_back("ATTTGG");
    // seqvec2.push_back("AGGTTT");

    // VDJRecombinationGenes genes("VA", alvec1, seqvec1, "JA", alvec2, seqvec2);

    // ClonotypeBuilder cl_builder;
    // // CCCG.AC.GGTTT
    // cl_builder.setSequence("CCCGACGGTTT")
    //         .setNucleotideSeq()
    //         .setRecombination(VJ_RECOMB)
    //         .addVarAlignment(1, 1, 1, 4)
    //         .addVarAlignment(3, 4, 3, 4)
    //         .addJoiAlignment(1, 2, 6, 5)
    //         .addJoiAlignment(2, 2, 9, 3)
    //         .addJoiAlignment(3, 2, 7, 5);
    // Clonotype clonotype = cl_builder.buildClonotype();
    // vector<Clonotype> vec;
    // vec.push_back(clonotype);
    // Cloneset cloneset(vec);

    // YMIR_ASSERT(writer.write(TEST_DATA_FOLDER + "../out.txt", cloneset, genes))

    // TODO: write a parser to test writer's output

YMIR_TEST_END


YMIR_TEST_START(test_parser_vj)

    // RepertoireParser parser;

    // bool V_err, J_err;
    // VDJRecombinationGenes vdj_genes("Vgene", TEST_DATA_FOLDER + "vgene.real.txt"
    //         , "Jgene", TEST_DATA_FOLDER + "jgene.real.txt", &V_err, &J_err);
    // YMIR_ASSERT(V_err)
    // YMIR_ASSERT(J_err)

    // Cloneset cr;
    // YMIR_ASSERT(parser.parse(TEST_DATA_FOLDER + "ymir.alpha.txt",
    //                          &cr,
    //                          vdj_genes,
    //                          NUCLEOTIDE,
    //                          VJ_RECOMB,
    //                          RepertoireParser::AlignmentColumnOptions()
    //                                  .setV(RepertoireParser::USE_PROVIDED)
    //                                  .setJ(RepertoireParser::USE_PROVIDED),
    //                          NaiveNucleotideAligner()))

    // YMIR_ASSERT(cr.size() == 30)
    // YMIR_ASSERT(cr[0].sequence() == "TGTGCAGCAAGTACCCCCTTAAGCTGGTGGTACTAGCTATGGAAAGCTGACATTT")
    // YMIR_ASSERT(cr[0].recombination() == VJ_RECOMB)
    // YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(0)].allele == "TRAV13-1")
    // YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(1)].allele == "TRAV13-2")
    // YMIR_ASSERT(cr[0].nVar() == 2)
    // YMIR_ASSERT(vdj_genes.J()[cr[0].getJoi(0)].allele == "TRAJ52")
    // YMIR_ASSERT(cr[0].nJoi() == 1)

    // YMIR_ASSERT(cr[2].sequence() == "TGTGCAACTCTTAGCAGGGATGAACACAGGCTTTCAGAAACTTGTATTT")
    // YMIR_ASSERT(cr[2].recombination() != VDJ_RECOMB)
    // YMIR_ASSERT(vdj_genes.V()[cr[2].getVar(0)].allele == "TRAV12-3")
    // YMIR_ASSERT(cr[2].nVar() == 1)
    // YMIR_ASSERT(vdj_genes.J()[cr[2].getJoi(0)].allele == "TRAJ8")
    // YMIR_ASSERT(vdj_genes.J()[cr[2].getJoi(1)].allele == "TRAJ18")
    // YMIR_ASSERT(cr[2].nJoi() == 2)

YMIR_TEST_END


YMIR_TEST_START(test_parser_vj_stream)

    // RepertoireParser parser;

    // bool V_err, J_err;
    // VDJRecombinationGenes vdj_genes("Vgene", TEST_DATA_FOLDER + "vgene.real.txt"
    //         , "Jgene", TEST_DATA_FOLDER + "jgene.real.txt", &V_err, &J_err);
    // YMIR_ASSERT(V_err)
    // YMIR_ASSERT(J_err)

    // Cloneset cr;
    // YMIR_ASSERT(parser.stream(TEST_DATA_FOLDER + "ymir.alpha.txt",
    //                          vdj_genes,
    //                          NUCLEOTIDE,
    //                          VJ_RECOMB,
    //                          RepertoireParser::AlignmentColumnOptions()
    //                                  .setV(RepertoireParser::USE_PROVIDED)
    //                                  .setJ(RepertoireParser::USE_PROVIDED),
    //                          NaiveNucleotideAligner()))

    // parser.nextBlock(&cr, 2);
    // YMIR_ASSERT(cr.size() == 2)
    // YMIR_ASSERT(cr[0].sequence() == "TGTGCAGCAAGTACCCCCTTAAGCTGGTGGTACTAGCTATGGAAAGCTGACATTT")
    // YMIR_ASSERT(cr[0].recombination() == VJ_RECOMB)
    // YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(0)].allele == "TRAV13-1")
    // YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(1)].allele == "TRAV13-2")
    // YMIR_ASSERT(cr[0].nVar() == 2)
    // YMIR_ASSERT(vdj_genes.J()[cr[0].getJoi(0)].allele == "TRAJ52")
    // YMIR_ASSERT(cr[0].nJoi() == 1)

    // parser.nextBlock(&cr, 5);
    // YMIR_ASSERT(cr.size() == 5)
    // YMIR_ASSERT(cr[0].sequence() == "TGTGCAACTCTTAGCAGGGATGAACACAGGCTTTCAGAAACTTGTATTT")
    // YMIR_ASSERT(cr[0].recombination() != VDJ_RECOMB)
    // YMIR_ASSERT(vdj_genes.V()[cr[2].getVar(0)].allele == "TRAV12-3")
    // YMIR_ASSERT(cr[0].nVar() == 1)
    // YMIR_ASSERT(vdj_genes.J()[cr[2].getJoi(0)].allele == "TRAJ8")
    // YMIR_ASSERT(vdj_genes.J()[cr[2].getJoi(1)].allele == "TRAJ18")
    // YMIR_ASSERT(cr[0].nJoi() == 2)

YMIR_TEST_END


YMIR_TEST_START(test_parser_vdj_with_d_alignment)

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

    // RepertoireParser parser;

    // Cloneset cr;
    // YMIR_ASSERT(parser.parse(TEST_DATA_FOLDER + "ymir.beta.txt",
    //                          &cr,
    //                          genes,
    //                          NUCLEOTIDE,
    //                          VDJ_RECOMB,
    //                          RepertoireParser::AlignmentColumnOptions()
    //                                  .setV(RepertoireParser::USE_PROVIDED)
    //                                  .setJ(RepertoireParser::USE_PROVIDED)
    //                                  .setD(RepertoireParser::OVERWRITE)))

    // YMIR_ASSERT2(cr.size(), 1)
    // YMIR_ASSERT2(cr[0].sequence(), "CCCGACGGTTT")
    // YMIR_ASSERT2(genes.V()[cr[0].getJoi(0)].allele, "Vseg1")
    // YMIR_ASSERT2(genes.J()[cr[0].getJoi(0)].allele, "Jseg1")
    // YMIR_ASSERT2( (int) cr[0].nDiv(), 3)
    // YMIR_ASSERT2( (int) cr[0].numDivAlignments(0), 1)
    // YMIR_ASSERT2( (int) cr[0].numDivAlignments(1), 2)
    // YMIR_ASSERT2( (int) cr[0].numDivAlignments(2), 3)
    // this is if default gene len is equal to 2
//    YMIR_ASSERT2( (int) cr[0].nDivAlignments(0), 3)
//    YMIR_ASSERT2( (int) cr[0].nDivAlignments(1), 4)
//    YMIR_ASSERT2( (int) cr[0].nDivAlignments(2), 5)

YMIR_TEST_END


YMIR_TEST_START(test_clorep)

    // RepertoireParser parser;

    // bool V_err, J_err;
    // VDJRecombinationGenes vdj_genes("Vgene", TEST_DATA_FOLDER + "vgene.real.txt"
    //         , "Jgene", TEST_DATA_FOLDER + "jgene.real.txt", &V_err, &J_err);
    // YMIR_ASSERT(V_err)
    // YMIR_ASSERT(J_err)

    // Cloneset cr;
    // YMIR_ASSERT(parser.parse(TEST_DATA_FOLDER + "ymir.alpha.txt",
    //                          &cr,
    //                          vdj_genes,
    //                          NUCLEOTIDE,
    //                          VJ_RECOMB,
    //                          RepertoireParser::AlignmentColumnOptions()
    //                                  .setV(RepertoireParser::USE_PROVIDED)
    //                                  .setJ(RepertoireParser::USE_PROVIDED)))

    // YMIR_ASSERT(!has_end_codon(cr[3].sequence()))
    // YMIR_ASSERT(is_out_of_frame(cr[3].sequence()))
    // YMIR_ASSERT(has_end_codon(cr[24].sequence()))
    // YMIR_ASSERT(is_out_of_frame(cr[24].sequence()))

    // ClonesetView crv = cr.head(10);
    // YMIR_ASSERT2(crv.size(), 10)
    // YMIR_ASSERT(cr[0].sequence() == crv[0].sequence())

    // vector<size_t> inds = {1, 5, 10};
    // crv = cr.subvec(inds);
    // YMIR_ASSERT(crv[0].sequence() == cr[1].sequence())
    // YMIR_ASSERT(crv[1].sequence() == cr[5].sequence())
    // YMIR_ASSERT(crv[2].sequence() == cr[10].sequence())

    // inds.clear(); inds.push_back(1); inds.push_back(2);
    // ClonesetView crv2 = crv.subvec(inds);

    // YMIR_ASSERT(crv2[0].sequence() == cr[5].sequence())
    // YMIR_ASSERT(crv2[1].sequence() == cr[10].sequence())

YMIR_TEST_END


YMIR_TEST_START(test_markovchain_nuc_mono)

    vector<prob_t> probs = {.1, .2, .3, .4};
    InsertionModel m(MONO_NUCLEOTIDE, probs.begin());

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
    m = InsertionModel(MONO_NUCLEOTIDE, probs.begin());
    std::default_random_engine rg;
    YMIR_ASSERT2(m.generate(5, rg), "AAAAA");

    probs = {0, 0, 0, 1};
    m = InsertionModel(MONO_NUCLEOTIDE, probs.begin());
    YMIR_ASSERT2(m.generate(1, rg), "T");

    probs = {0, 1, 0, 0};
    m = InsertionModel(MONO_NUCLEOTIDE, probs.begin());
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

    InsertionModel mc(mat);
    string s = "ACGT";

    // .25 * .7 * .3 * .2 = .0105
    YMIR_ASSERT(mc.nucProbability("") == 1);
    YMIR_ASSERT2(mc.nucProbability(s) - .0105, 0);
    YMIR_ASSERT(mc.nucProbability(s.begin(), 4) == .0105);

    // .4 * .7 * .3 * .2 = .0126
    YMIR_ASSERT2(mc.nucProbability(s.begin(), 4, 'C'), .0168);

    mc.updateProbabilities(vec.begin());

    YMIR_ASSERT2(mc.nucProbability(s) - .0105, 0);
    YMIR_ASSERT2(mc.nucProbability(s.begin(), 4, 'C'), .0168);

    // .25 * .25 * .1
    YMIR_ASSERT2(mc.nucProbability(s.rbegin(), 3, '_'), .00625);

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

    ClonotypeBuilder cl_builder;
    // CCCG.AC.GGTTT
    // poses:
    // vs: 0-1-2-3-4-5
    // js: 7-8-9-10-11-12
    cl_builder.setSequence("CCCGACGGTTT")
            .setNucleotideSeq()
            .addVarAlignment(1, 1, 1, 4)
            .addVarAlignment(3, 1, 1, 5)
            .addJoiAlignment(1, 3, 8, 4)
            .addJoiAlignment(2, 2, 9, 3)
            .addJoiAlignment(3, 2, 7, 5);
    cl_builder.setRecombination(VJ_RECOMB);
    Clonotype clonotype = cl_builder.buildClonotype();

//    cout << "here" << endl;
    MAAG maag = maag_builder.build(clonotype, SAVE_METADATA);

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

    YMIR_ASSERT2(maag.event_index(3, 2, 0, 0), mvec.event_index(VJ_JOI_DEL, 2, 1))
    YMIR_ASSERT2(maag.event_index(3, 2, 1, 0), mvec.event_index(VJ_JOI_DEL, 2, 2))
    YMIR_ASSERT2(maag.event_index(3, 2, 2, 0), mvec.event_index(VJ_JOI_DEL, 2, 3))
    YMIR_ASSERT2(maag.event_index(3, 2, 3, 0), mvec.event_index(VJ_JOI_DEL, 2, 4))
    YMIR_ASSERT2(maag.event_index(3, 2, 4, 0), mvec.event_index(VJ_JOI_DEL, 2, 5))
    YMIR_ASSERT2(maag.event_index(3, 2, 5, 0), mvec.event_index(VJ_JOI_DEL, 2, 6))

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

    ClonotypeBuilder cl_builder;
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
            .setNucleotideSeq()
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
    Clonotype clonotype = cl_builder.buildClonotype();

    MAAG maag = maag_builder.build(clonotype, SAVE_METADATA);

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
    YMIR_ASSERT2(maag.fullProbability(0, 0, 0), maag_builder.build(clonotype, NO_METADATA).fullProbability(0, 0, 0))  // error is here
    YMIR_ASSERT2(maag.fullProbability(1, 1, 1), maag_builder.build(clonotype, NO_METADATA).fullProbability(1, 1, 1))
    YMIR_ASSERT2(maag.fullProbability(0, 2, 2), maag_builder.build(clonotype, NO_METADATA).fullProbability(0, 2, 2))  // error is here

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

    ClonotypeBuilder cl_builder;
    // CCCG.AC.GGTTT
    // poses:
    // vs: 0-1-2-3-4-5
    // js: 7-8-9-10-11
    cl_builder.setSequence("CCCGACGGTTT")
            .setNucleotideSeq()
            .addVarAlignment(1, 1, 1, 4)
            .addVarAlignment(3, 1, 1, 5)
            .addJoiAlignment(1, 2, 8, 4)
            .addJoiAlignment(2, 2, 9, 3)
            .addJoiAlignment(3, 2, 7, 5);
    cl_builder.setRecombination(VJ_RECOMB);
    Clonotype clonotype = cl_builder.buildClonotype();

    MAAG maag = maag_builder.build(clonotype, SAVE_METADATA);

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

    ClonotypeBuilder cl_builder;
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
            .setNucleotideSeq()
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
    Clonotype clonotype = cl_builder.buildClonotype();

    MAAG maag = maag_builder.build(clonotype, SAVE_METADATA);

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


YMIR_TEST_START(test_model_vj_file)

    ProbabilisticAssemblingModel model1(TEST_DATA_FOLDER + "randomfile");
    YMIR_ASSERT(!model1.status())

    ModelParameterVector mvec = make_test_events_vj();

    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vj_model/");
    YMIR_ASSERT(model.status())

    YMIR_ASSERT(mvec == model.event_probabilities())

YMIR_TEST_END


YMIR_TEST_START(test_model_vdj_file)

    ModelParameterVector mvec = make_test_events_vdj();

    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vdj_model/");
    YMIR_ASSERT(model.status())

    YMIR_ASSERT(mvec == model.event_probabilities())

YMIR_TEST_END


YMIR_TEST_START(test_model_vj_save_load)

    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vj_model/");
    YMIR_ASSERT(model.status())

    YMIR_ASSERT(model.save(TEST_DATA_FOLDER + "test_vj_model2/"))
    ProbabilisticAssemblingModel model2(TEST_DATA_FOLDER + "test_vj_model2/");
    YMIR_ASSERT(model2.status())

    YMIR_ASSERT(model.event_probabilities() == model2.event_probabilities())

YMIR_TEST_END


YMIR_TEST_START(test_model_vdj_save_load)

    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vdj_model/");
    YMIR_ASSERT(model.status())

    YMIR_ASSERT(model.save(TEST_DATA_FOLDER + "test_vdj_model2/"))
    ProbabilisticAssemblingModel model2(TEST_DATA_FOLDER + "test_vdj_model2/");
    YMIR_ASSERT(model2.status())

    YMIR_ASSERT(model.event_probabilities() == model2.event_probabilities())

YMIR_TEST_END


YMIR_TEST_START(test_model_gene_usage)

    ProbabilisticAssemblingModel model_vj(TEST_DATA_FOLDER + "test_vj_model/");
    YMIR_ASSERT(model_vj.status())

    ProbabilisticAssemblingModel model_vdj(TEST_DATA_FOLDER + "test_vdj_model/");
    YMIR_ASSERT(model_vdj.status())

YMIR_TEST_END


YMIR_TEST_START(test_model_vj_maag)

    // TGTGCTCTTGGGGAACTTTCGGAGTGGCTCTAGCAACACAGGCAAACTAATCTTT	CALGELSEW~SSNTGKLIF	TRAV13-2		TRAJ37	1|1|5		0|25|31

    ModelParameterVector mvec = make_test_events_vj();

    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vj_model/");
    YMIR_ASSERT(model.status())

    // RepertoireParser parser;

    // bool V_err, J_err;
    // VDJRecombinationGenes vj_genes("Vgene", TEST_DATA_FOLDER + "vgene.real.txt",
    //                                 "Jgene", TEST_DATA_FOLDER + "jgene.real.txt",
    //                                 &V_err, &J_err);
    // YMIR_ASSERT(V_err)
    // YMIR_ASSERT(J_err)

    // Cloneset cloneset;
    // YMIR_ASSERT(parser.parse(TEST_DATA_FOLDER + "ymir.alpha2.txt", &cloneset, vj_genes, NUCLEOTIDE, VJ_RECOMB))

    // MAAG maag = model.buildGraphs(cloneset, SAVE_METADATA, NUCLEOTIDE)[1];

    // YMIR_ASSERT2(maag.event_index(0, 0, 0, 0), mvec.event_index(VJ_VAR_JOI_GEN, 0, 0, 0))
    // YMIR_ASSERT2(maag.event_index(0, 0, 0, 1), mvec.event_index(VJ_VAR_JOI_GEN, 0, 0, 1))
    // YMIR_ASSERT2(maag.event_index(0, 0, 1, 0), mvec.event_index(VJ_VAR_JOI_GEN, 0, 2, 0))
    // YMIR_ASSERT2(maag.event_index(0, 0, 1, 2), mvec.event_index(VJ_VAR_JOI_GEN, 0, 2, 2))

    // YMIR_ASSERT2(maag.event_index(1, 0, 0, 0), mvec.event_index(VJ_VAR_DEL, 0, 4));
    // YMIR_ASSERT2(maag.event_index(1, 0, 0, 1), mvec.event_index(VJ_VAR_DEL, 0, 3));
    // YMIR_ASSERT2(maag.event_index(1, 0, 0, 2), mvec.event_index(VJ_VAR_DEL, 0, 2));
    // YMIR_ASSERT2(maag.event_index(1, 0, 0, 3), mvec.event_index(VJ_VAR_DEL, 0, 1));
    // YMIR_ASSERT2(maag.event_index(1, 0, 0, 4), mvec.event_index(VJ_VAR_DEL, 0, 0));
    // YMIR_ASSERT2(maag.event_index(1, 0, 0, 5), 0)

    // YMIR_ASSERT2(maag.event_index(1, 1, 0, 0), mvec.event_index(VJ_VAR_DEL, 2, 6))
    // YMIR_ASSERT2(maag.event_index(1, 1, 0, 1), mvec.event_index(VJ_VAR_DEL, 2, 5))
    // YMIR_ASSERT2(maag.event_index(1, 1, 0, 2), mvec.event_index(VJ_VAR_DEL, 2, 4))
    // YMIR_ASSERT2(maag.event_index(1, 1, 0, 3), mvec.event_index(VJ_VAR_DEL, 2, 3))
    // YMIR_ASSERT2(maag.event_index(1, 1, 0, 4), mvec.event_index(VJ_VAR_DEL, 2, 2))
    // YMIR_ASSERT2(maag.event_index(1, 1, 0, 5), mvec.event_index(VJ_VAR_DEL, 2, 1))

    // YMIR_ASSERT2(maag.event_index(2, 0, 0, 0), mvec.event_index(VJ_VAR_JOI_INS_LEN, 0, maag.position(6) - maag.position(0) - 1))
    // YMIR_ASSERT2(maag.event_index(2, 0, 0, 1), 0)

    // YMIR_ASSERT2(maag.event_index(3, 0, 0, 0), 0)
    // YMIR_ASSERT2(maag.event_index(3, 0, 1, 0), mvec.event_index(VJ_JOI_DEL, 0, 2))
    // YMIR_ASSERT2(maag.event_index(3, 0, 2, 0), mvec.event_index(VJ_JOI_DEL, 0, 3))
    // YMIR_ASSERT2(maag.event_index(3, 0, 3, 0), mvec.event_index(VJ_JOI_DEL, 0, 4))
    // YMIR_ASSERT2(maag.event_index(3, 0, 4, 0), mvec.event_index(VJ_JOI_DEL, 0, 5))
    // YMIR_ASSERT2(maag.event_index(3, 0, 5, 0), mvec.event_index(VJ_JOI_DEL, 0, 6))

    // YMIR_ASSERT2(maag.event_index(3, 2, 0, 0), mvec.event_index(VJ_JOI_DEL, 2, 1))
    // YMIR_ASSERT2(maag.event_index(3, 2, 1, 0), mvec.event_index(VJ_JOI_DEL, 2, 2))
    // YMIR_ASSERT2(maag.event_index(3, 2, 2, 0), mvec.event_index(VJ_JOI_DEL, 2, 3))
    // YMIR_ASSERT2(maag.event_index(3, 2, 3, 0), mvec.event_index(VJ_JOI_DEL, 2, 4))
    // YMIR_ASSERT2(maag.event_index(3, 2, 4, 0), mvec.event_index(VJ_JOI_DEL, 2, 5))
    // YMIR_ASSERT2(maag.event_index(3, 2, 5, 0), mvec.event_index(VJ_JOI_DEL, 2, 6))

    // YMIR_ASSERT2(maag.event_probability(0, 0, 0, 0), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 0, 0))
    // YMIR_ASSERT2(maag.event_probability(0, 0, 0, 1), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 0, 1))
    // YMIR_ASSERT2(maag.event_probability(0, 0, 1, 0), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 2, 0))
    // YMIR_ASSERT2(maag.event_probability(0, 0, 1, 2), mvec.event_prob(VJ_VAR_JOI_GEN, 0, 2, 2))

    // YMIR_ASSERT2(maag.event_probability(1, 0, 0, 0), mvec.event_prob(VJ_VAR_DEL, 0, 4));
    // YMIR_ASSERT2(maag.event_probability(1, 0, 0, 1), mvec.event_prob(VJ_VAR_DEL, 0, 3));
    // YMIR_ASSERT2(maag.event_probability(1, 0, 0, 2), mvec.event_prob(VJ_VAR_DEL, 0, 2));
    // YMIR_ASSERT2(maag.event_probability(1, 0, 0, 3), mvec.event_prob(VJ_VAR_DEL, 0, 1));
    // YMIR_ASSERT2(maag.event_probability(1, 0, 0, 4), mvec.event_prob(VJ_VAR_DEL, 0, 0));
    // YMIR_ASSERT2(maag.event_probability(1, 0, 0, 5), 0)

    // YMIR_ASSERT2(maag.event_probability(1, 1, 0, 0), mvec.event_prob(VJ_VAR_DEL, 2, 6))
    // YMIR_ASSERT2(maag.event_probability(1, 1, 0, 1), mvec.event_prob(VJ_VAR_DEL, 2, 5))
    // YMIR_ASSERT2(maag.event_probability(1, 1, 0, 2), mvec.event_prob(VJ_VAR_DEL, 2, 4))
    // YMIR_ASSERT2(maag.event_probability(1, 1, 0, 3), mvec.event_prob(VJ_VAR_DEL, 2, 3))
    // YMIR_ASSERT2(maag.event_probability(1, 1, 0, 4), mvec.event_prob(VJ_VAR_DEL, 2, 2))
    // YMIR_ASSERT2(maag.event_probability(1, 1, 0, 5), mvec.event_prob(VJ_VAR_DEL, 2, 1))

    // YMIR_ASSERT2(maag.event_probability(2, 0, 0, 1), 0)

    // YMIR_ASSERT2(maag.event_probability(3, 0, 0, 0), 0)
    // YMIR_ASSERT2(maag.event_probability(3, 0, 1, 0), mvec.event_prob(VJ_JOI_DEL, 0, 2))
    // YMIR_ASSERT2(maag.event_probability(3, 0, 2, 0), mvec.event_prob(VJ_JOI_DEL, 0, 3))
    // YMIR_ASSERT2(maag.event_probability(3, 0, 3, 0), mvec.event_prob(VJ_JOI_DEL, 0, 4))
    // YMIR_ASSERT2(maag.event_probability(3, 0, 4, 0), mvec.event_prob(VJ_JOI_DEL, 0, 5))
    // YMIR_ASSERT2(maag.event_probability(3, 0, 5, 0), mvec.event_prob(VJ_JOI_DEL, 0, 6))

    // YMIR_ASSERT2(maag.event_probability(3, 2, 0, 0), mvec.event_prob(VJ_JOI_DEL, 2, 1))
    // YMIR_ASSERT2(maag.event_probability(3, 2, 1, 0), mvec.event_prob(VJ_JOI_DEL, 2, 2))
    // YMIR_ASSERT2(maag.event_probability(3, 2, 2, 0), mvec.event_prob(VJ_JOI_DEL, 2, 3))
    // YMIR_ASSERT2(maag.event_probability(3, 2, 3, 0), mvec.event_prob(VJ_JOI_DEL, 2, 4))
    // YMIR_ASSERT2(maag.event_probability(3, 2, 4, 0), mvec.event_prob(VJ_JOI_DEL, 2, 5))
    // YMIR_ASSERT2(maag.event_probability(3, 2, 5, 0), mvec.event_prob(VJ_JOI_DEL, 2, 6))

YMIR_TEST_END


YMIR_TEST_START(test_model_vdj_maag)

    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vdj_model/");
    YMIR_ASSERT(model.status())

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_maag_forward_backward_vj)

    ModelParameterVector mvec = make_test_events_vj();

    vector<string> alvec1;
    vector<string> seqvec1;
    alvec1.push_back("Vseg1");
    alvec1.push_back("Vseg2");
    alvec1.push_back("Vseg3");
    seqvec1.push_back("CCCA");
    seqvec1.push_back("GGG");
    seqvec1.push_back("CCCAGG");

    vector<string> alvec2;
    vector<string> seqvec2;
    alvec2.push_back("Jseg1");
    alvec2.push_back("Jseg2");
    alvec2.push_back("Jseg3");
    seqvec2.push_back("CCGTTT");
    seqvec2.push_back("CATT");
    seqvec2.push_back("AGGTTT");

    VDJRecombinationGenes genes("VA", alvec1, seqvec1, "JA", alvec2, seqvec2);

    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 0)] = 1;
    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 1)] = 0;
    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 2)] = 0;
    mvec[mvec.event_index(VJ_VAR_JOI_INS_NUC, 0, 3)] = 0;
    MAAGBuilder maag_builder(mvec, genes);

    ClonotypeBuilder cl_builder;
    // CCCG.AC.GGTTT
    // poses:
    // vs: 0-1-2-3-4-5
    // js: 7-8-9-10-11
    cl_builder.setSequence("CCCAAAAAAATT")
            .setNucleotideSeq()
            .addVarAlignment(1, 1, 1, 4)
//            .addVarAlignment(3, 5)
//            .addJoiAlignment(1, 8)
//            .addJoiAlignment(2, 9)
            .addJoiAlignment(2, 2, 10, 3);
//            .addJoiAlignment(3, 7);
    cl_builder.setRecombination(VJ_RECOMB);
    Clonotype clonotype = cl_builder.buildClonotype();

    MAAG maag = maag_builder.build(clonotype, SAVE_METADATA);

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
    Clonotype clonotype = cl_builder.buildClonotype();

    MAAG maag = maag_builder.build(clonotype, SAVE_METADATA);

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


struct TestInfo {
    string test_name;
    vector<string> failed_cases;

    TestInfo(const string& name, const vector<string> vec) : test_name(name), failed_cases(vec) {}
};


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
    YMIR_TEST(test_basic())

    // Tests for ModelParameterVector
    YMIR_TEST(test_model_param_vec_vj())
    YMIR_TEST(test_model_param_vec_vdj())

    // Tests for gene segments classes
    YMIR_TEST(test_genesegmentalphabet())
    YMIR_TEST(test_vdjgenes2())
    YMIR_TEST(test_vdjgenes3())
    YMIR_TEST(test_genesegmentalphabet_read())
    YMIR_TEST(test_vdjgenes_read())

    // Tests for clone, clone alignment and clone builder classes.
    YMIR_TEST(test_clone())
    YMIR_TEST(test_clonebuilder_clonealign())

    // Tests for NoGapAlignmentVector and GappedAlignmentVector
    YMIR_TEST(test_nogap_alignment_vector_no_err())
    YMIR_TEST(test_nogap_alignment_vector_errors())
    YMIR_TEST(test_gapped_alignment_vector())

    // Tests for VDJAlignment and VDJAlignmentBuilder
    YMIR_TEST(test_vdj_alignment_simple_vj())
    YMIR_TEST(test_vdj_alignment_simple_vdj())
    YMIR_TEST(test_vdj_alignment_vector_vj())
    YMIR_TEST(test_vdj_alignment_vector_vdj())

    // Tests for sequences aligners.
    YMIR_TEST(test_naive_cdr3_nuc_aligner())
    YMIR_TEST(test_cdr3_nuc_aligner())
    YMIR_TEST(test_cdr3_aa_aligner())
    YMIR_TEST(test_sw_aligner())
    YMIR_TEST(test_swng_aligner())

    // Error corrector test
    YMIR_TEST(test_errcorr_aligner())

    // Test for the repertoire parser and writer
    // YMIR_TEST(test_parser_vj())
    // YMIR_TEST(test_parser_vj_stream())
    // YMIR_TEST(test_parser_vdj_with_d_alignment())
//    YMIR_TEST(test_ymir_vdj_wo_d_alignment())
    YMIR_TEST(test_writer())

    // Tests for clonal repertoires and clonal repertoire views.
    YMIR_TEST(test_clorep())

    // Tests for markov chain.
    YMIR_TEST(test_markovchain_nuc_mono())
    YMIR_TEST(test_markovchain_nuc_di())
    YMIR_TEST(test_markovchain_aa())

    // Test for Multi-Matrix Chains
    YMIR_TEST(test_mmc())

    // Tests for MAAG / MAAG builder
    YMIR_TEST(test_maag_vj())
    YMIR_TEST(test_maag_vdj())
    YMIR_TEST(test_maag_builder_replace_vj())
    YMIR_TEST(test_maag_builder_replace_vdj())

    // Tests for probabilistic assembling model (PAM) reading / writing files.
    YMIR_TEST(test_model_vj_file())
    YMIR_TEST(test_model_vdj_file())
    YMIR_TEST(test_model_vj_save_load())
    YMIR_TEST(test_model_vdj_save_load())
    YMIR_TEST(test_model_gene_usage())
//    YMIR_TEST(test_model_vj_maag())
//    YMIR_TEST(test_model_vdj_maag())

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