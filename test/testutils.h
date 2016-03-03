//
// Created by Vadim N. on 20/02/2016.
//

#ifndef YMIR_TESTUTILS_H_H
#define YMIR_TESTUTILS_H_H


#include <iostream>
#include <list>
#include <sstream>

#include "Inference"


namespace ymir {

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

    struct TestInfo {
        string test_name;
        vector<string> failed_cases;

        TestInfo(const string& name, const vector<string> vec) : test_name(name), failed_cases(vec) {}
    };

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

}

#endif //YMIR_TESTUTILS_H_H
