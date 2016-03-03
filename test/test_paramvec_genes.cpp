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
    
    // Tests for ModelParameterVector
    YMIR_TEST(test_model_param_vec_vj())
    YMIR_TEST(test_model_param_vec_vdj())

    // Tests for gene segments classes
    YMIR_TEST(test_genesegmentalphabet())
    YMIR_TEST(test_vdjgenes2())
    YMIR_TEST(test_vdjgenes3())
    YMIR_TEST(test_genesegmentalphabet_read())
    YMIR_TEST(test_vdjgenes_read())

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