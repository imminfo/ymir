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


YMIR_TEST_START(test_codon_alignment_vector)
    YMIR_ASSERT(false)
YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_nuc_simple_vj)

    VDJAlignmentBuilder builder;

    builder.addVarAlignment(11, 1, 2, 3)
           .addVarAlignment(12, 4, 5, 6)
           .addJoiAlignment(31, 7, 8, 9)
           .addJoiAlignment(33, 10, 11, 12)
           .addJoiAlignment(35, 13, 14, 15);

    VDJAlignment algn = builder.buildAlignment();

    YMIR_ASSERT2(algn.nVar(), 2)
    YMIR_ASSERT2(algn.nJoi(), 3)
    YMIR_ASSERT2(algn.nDiv(), 0)

    YMIR_ASSERT2(algn.getVar(0), 11)
    YMIR_ASSERT2(algn.getVar(1), 12)
    YMIR_ASSERT2(algn.getJoi(0), 31)
    YMIR_ASSERT2(algn.getJoi(1), 33)
    YMIR_ASSERT2(algn.getJoi(2), 35)

    YMIR_ASSERT2(algn.getVarGeneStart(0), 1)
    YMIR_ASSERT2(algn.getVarGeneEnd(0), 3)
    YMIR_ASSERT2(algn.getVarSeqStart(0), 2)
    YMIR_ASSERT2(algn.getVarSeqEnd(0), 4)
    YMIR_ASSERT2(algn.getVarLen(0), 3)
    YMIR_ASSERT2(algn.getVarGeneStart(1), 4)
    YMIR_ASSERT2(algn.getVarGeneEnd(1), 9)
    YMIR_ASSERT2(algn.getVarSeqStart(1), 5)
    YMIR_ASSERT2(algn.getVarSeqEnd(1), 10)
    YMIR_ASSERT2(algn.getVarLen(1), 6)

    YMIR_ASSERT2(algn.getJoiGeneStart(0), 7)
    YMIR_ASSERT2(algn.getJoiGeneEnd(0), 15)
    YMIR_ASSERT2(algn.getJoiSeqStart(0), 8)
    YMIR_ASSERT2(algn.getJoiSeqEnd(0), 16)
    YMIR_ASSERT2(algn.getJoiLen(0), 9)
    YMIR_ASSERT2(algn.getJoiGeneStart(2), 13)
    YMIR_ASSERT2(algn.getJoiGeneEnd(2), 27)
    YMIR_ASSERT2(algn.getJoiSeqStart(2), 14)
    YMIR_ASSERT2(algn.getJoiSeqEnd(2), 28)
    YMIR_ASSERT2(algn.getJoiLen(2), 15)


    builder.addVarAlignment(41, 41, 42, 43)
           .addVarAlignment(42, 44, 45, 46)
           .addVarAlignment(49, 57, 58, 59)
           .addJoiAlignment(53, 1, 2, 3)
           .addJoiAlignment(55, 4, 5, 6);

    algn = builder.buildAlignment();

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


YMIR_TEST_START(test_vdj_alignment_nuc_simple_vdj)

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

    VDJAlignment algn = builder.buildAlignment();

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
    YMIR_ASSERT2(algn.getDivGeneEnd(0, 0), 41)
    YMIR_ASSERT2(algn.getDivSeqStart(0, 0), 21)
    YMIR_ASSERT2(algn.getDivSeqEnd(0, 0), 42)
    YMIR_ASSERT2(algn.getDivLen(0, 0), 22)
    YMIR_ASSERT2(algn.getDivGeneStart(0, 1), 23)
    YMIR_ASSERT2(algn.getDivGeneEnd(0, 1), 47)
    YMIR_ASSERT2(algn.getDivSeqStart(0, 1), 24)
    YMIR_ASSERT2(algn.getDivSeqEnd(0, 1), 48)
    YMIR_ASSERT2(algn.getDivLen(0, 1), 25)

    YMIR_ASSERT2(algn.getDivGeneStart(1, 0), 30)
    YMIR_ASSERT2(algn.getDivGeneEnd(1, 0), 61)
    YMIR_ASSERT2(algn.getDivSeqStart(1, 0), 31)
    YMIR_ASSERT2(algn.getDivSeqEnd(1, 0), 62)
    YMIR_ASSERT2(algn.getDivLen(1, 0), 32)
    YMIR_ASSERT2(algn.getDivGeneStart(1, 2), 36)
    YMIR_ASSERT2(algn.getDivGeneEnd(1, 2), 73)
    YMIR_ASSERT2(algn.getDivSeqStart(1, 2), 37)
    YMIR_ASSERT2(algn.getDivSeqEnd(1, 2), 74)
    YMIR_ASSERT2(algn.getDivLen(1, 2), 38)

YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_nuc_vector_vj)

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

    VDJAlignment algn = builder.buildAlignment();

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


YMIR_TEST_START(test_vdj_alignment_nuc_vector_vdj)

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

    VDJAlignment algn = builder.buildAlignment();

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
    YMIR_ASSERT2(algn.getDiv(0), 43)
    YMIR_ASSERT2(algn.getDiv(1), 45)

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

     YMIR_ASSERT2(algn.getDivGeneStart(0, 0), 13)
     YMIR_ASSERT2(algn.getDivSeqStart(0, 0), 14)
     YMIR_ASSERT2(algn.getDivLen(0, 0), 6)
     YMIR_ASSERT(!algn.isDivMismatch(0, 0, 1))
     YMIR_ASSERT(algn.isDivMismatch(0, 0, 2))
     YMIR_ASSERT(algn.isDivMismatch(0, 0, 3))
     YMIR_ASSERT(!algn.isDivMismatch(0, 0, 4))
     YMIR_ASSERT(algn.isDivMismatch(0, 0, 5))
     YMIR_ASSERT(algn.isDivMismatch(0, 0, 6))

    YMIR_ASSERT2(algn.numDivMismatches(0, 0, 1, 6), 4)
    YMIR_ASSERT2(algn.numDivMismatches(0, 0, 3, 6), 3)
    YMIR_ASSERT2(algn.numDivMismatches(0, 0, 2, 5), 3)

     YMIR_ASSERT2(algn.getDivGeneStart(0, 1), 15)
     YMIR_ASSERT2(algn.getDivSeqStart(0, 1), 16)
     YMIR_ASSERT2(algn.getDivLen(0, 1), 3)
     YMIR_ASSERT(!algn.isDivMismatch(0, 1, 1))
     YMIR_ASSERT(algn.isDivMismatch(0, 1, 2))
     YMIR_ASSERT(!algn.isDivMismatch(0, 1, 3))

     YMIR_ASSERT2(algn.getDivGeneStart(1, 0), 17)
     YMIR_ASSERT2(algn.getDivSeqStart(1, 0), 18)
     YMIR_ASSERT2(algn.getDivLen(1, 0), 7)
     YMIR_ASSERT(!algn.isDivMismatch(1, 0, 1))
     YMIR_ASSERT(!algn.isDivMismatch(1, 0, 2))
     YMIR_ASSERT(algn.isDivMismatch(1, 0, 3))
     YMIR_ASSERT(algn.isDivMismatch(1, 0, 4))
     YMIR_ASSERT(!algn.isDivMismatch(1, 0, 5))
     YMIR_ASSERT(algn.isDivMismatch(1, 0, 6))
     YMIR_ASSERT(algn.isDivMismatch(1, 0, 7))

YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_aa_simple_vj())

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_aa_simple_vdj())

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_aa_vector_vj())

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_vdj_alignment_aa_vector_vdj())

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_clone)

    VDJAlignmentBuilder builder;

    builder.addVarAlignment(1, 8, 9, 5);
    builder.addVarAlignment(2, 11, 12, 3);
    builder.addVarAlignment(3, 14, 15, 6);

    builder.addJoiAlignment(4, 17, 18, 7);
    builder.addJoiAlignment(5, 20, 21, 3);

    builder.addDivAlignment(6, 31, 32, 4);
    builder.addDivAlignment(6, 41, 42, 7);

    builder.addDivAlignment(7, 51, 52, 6);
    builder.addDivAlignment(7, 60, 61, 4);

    Clonotype c("cloneseq", NUCLEOTIDE, VDJ_RECOMB, builder.buildAlignment());
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

    YMIR_ASSERT2(c.getVarGeneStart(0), 8)
    YMIR_ASSERT2(c.getVarSeqEnd(0), 13)
    YMIR_ASSERT2(c.getVarLen(0), 5)

    YMIR_ASSERT2(c.getVarGeneEnd(2), 19)
    YMIR_ASSERT2(c.getVarSeqStart(2), 15)
    YMIR_ASSERT2(c.getVarLen(2), 6)

    YMIR_ASSERT2(c.getJoiGeneStart(0), 17)
    YMIR_ASSERT2(c.getJoiSeqEnd(0), 24)
    YMIR_ASSERT2(c.getJoiLen(0), 7)

    YMIR_ASSERT2(c.getJoiGeneEnd(1), 22)
    YMIR_ASSERT2(c.getJoiSeqStart(1), 21)
    YMIR_ASSERT2(c.getJoiLen(1), 3)

    YMIR_ASSERT2(c.getDivGeneStart(0, 0), 31)
    YMIR_ASSERT2(c.getDivGeneEnd(0, 0), 34)
    YMIR_ASSERT2(c.getDivSeqStart(0, 0), 32)
    YMIR_ASSERT2(c.getDivSeqEnd(0, 0), 35)

    YMIR_ASSERT2(c.getDivSeqStart(0, 1), 42)
    YMIR_ASSERT2(c.getDivSeqEnd(0, 1), 48)

    YMIR_ASSERT2(c.getDivGeneStart(1, 0), 51)
    YMIR_ASSERT2(c.getDivGeneEnd(1, 0), 56)
    YMIR_ASSERT2(c.getDivSeqStart(1, 0), 52)
    YMIR_ASSERT2(c.getDivSeqEnd(1, 0), 57)

    YMIR_ASSERT2(c.getDivGeneStart(1, 1), 60)
    YMIR_ASSERT2(c.getDivGeneEnd(1, 1), 63)
    YMIR_ASSERT2(c.getDivSeqStart(1, 1), 61)
    YMIR_ASSERT2(c.getDivSeqEnd(1, 1), 64)
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
    YMIR_ASSERT2(c.getVarGeneStart(0), 1)
    YMIR_ASSERT2(c.getVarSeqStart(0), 3)
    YMIR_ASSERT2(c.getVarLen(0), 15)

    YMIR_ASSERT2(c.getVar(1), 11)
    YMIR_ASSERT2(c.getVarGeneStart(1), 4)
    YMIR_ASSERT2(c.getVarSeqStart(1), 1)
    YMIR_ASSERT2(c.getVarLen(1), 25)

    YMIR_ASSERT2(c.getVar(2), 12)
    YMIR_ASSERT2(c.getVarGeneStart(2), 2)
    YMIR_ASSERT2(c.getVarSeqStart(2), 6)
    YMIR_ASSERT2(c.getVarLen(2), 35)

    YMIR_ASSERT2(c.getJoi(0), 20)
    YMIR_ASSERT2(c.getJoiGeneStart(0), 1)
    YMIR_ASSERT2(c.getJoiSeqStart(0), 21)
    YMIR_ASSERT2(c.getJoiLen(0), 10)

    YMIR_ASSERT2(c.getDiv(0), 31)
    YMIR_ASSERT2(c.getDivGeneStart(0, 0), 8)
    YMIR_ASSERT2(c.getDivSeqStart(0, 0), 11)
    YMIR_ASSERT2(c.getDivLen(0, 0), 2)
    YMIR_ASSERT2(c.getDivGeneStart(0, 1), 8)
    YMIR_ASSERT2(c.getDivSeqStart(0, 1), 13)
    YMIR_ASSERT2(c.getDivLen(0, 1), 2)

    YMIR_ASSERT2(c.getDiv(1), 30)
    YMIR_ASSERT2(c.getDivGeneStart(1, 0), 1)
    YMIR_ASSERT2(c.getDivSeqStart(1, 0), 3)
    YMIR_ASSERT2(c.getDivLen(1, 0), 2)
    YMIR_ASSERT2(c.getDivGeneStart(1, 1), 1)
    YMIR_ASSERT2(c.getDivSeqStart(1, 1), 5)
    YMIR_ASSERT2(c.getDivLen(1, 1), 2)
    YMIR_ASSERT2(c.getDiv(2), 32)
    YMIR_ASSERT2(c.getDivGeneStart(2, 0), 1)
    YMIR_ASSERT2(c.getDivSeqStart(2, 0), 5)
    YMIR_ASSERT2(c.getDivLen(2, 0), 2)

    YMIR_ASSERT2(c.nVar(), 3)
    YMIR_ASSERT2(c.nJoi(), 1)
    YMIR_ASSERT2(c.nDiv(), 3)
    YMIR_ASSERT2(c.numDivAlignments(0), 2)
    YMIR_ASSERT2(c.numDivAlignments(1), 2)
    YMIR_ASSERT2(c.numDivAlignments(2), 1)

YMIR_TEST_END


YMIR_TEST_START(test_clorep)

     NaiveNucParser parser;

     bool V_err, J_err;
     VDJRecombinationGenes vdj_genes("Vgene", TEST_DATA_FOLDER + "vgene.real.txt"
             , "Jgene", TEST_DATA_FOLDER + "jgene.real.txt", &V_err, &J_err);
     YMIR_ASSERT(V_err)
     YMIR_ASSERT(J_err)

     Cloneset cr;
     YMIR_ASSERT(parser.openAndParse(TEST_DATA_FOLDER + "ymir.alpha.txt",
                                     &cr,
                                     vdj_genes,
                                     NUCLEOTIDE,
                                     VJ_RECOMB,
                                     AlignmentColumnOptions()
                                      .setV(AlignmentColumnOptions::USE_PROVIDED)
                                      .setJ(AlignmentColumnOptions::USE_PROVIDED)))

    YMIR_ASSERT(!has_end_codon(cr[3].sequence()))
    YMIR_ASSERT(cr[3].isNoncoding())
    YMIR_ASSERT(is_out_of_frame(cr[3].sequence()))
    YMIR_ASSERT(cr[3].isOutOfFrame())
    YMIR_ASSERT(has_end_codon(cr[24].sequence()))
    YMIR_ASSERT(cr[24].isNoncoding())
    YMIR_ASSERT(is_out_of_frame(cr[24].sequence()))
    YMIR_ASSERT(cr[24].isOutOfFrame())

     ClonesetView crv = cr.head(10);
     YMIR_ASSERT2(crv.size(), 10)
     YMIR_ASSERT(cr[0].sequence() == crv[0].sequence())

     vector<size_t> inds = {1, 5, 10};
     crv = cr.subvec(inds);
     YMIR_ASSERT(crv[0].sequence() == cr[1].sequence())
     YMIR_ASSERT(crv[1].sequence() == cr[5].sequence())
     YMIR_ASSERT(crv[2].sequence() == cr[10].sequence())

     inds.clear(); inds.push_back(1); inds.push_back(2);
     ClonesetView crv2 = crv.subvec(inds);

     YMIR_ASSERT(crv2[0].sequence() == cr[5].sequence())
     YMIR_ASSERT(crv2[1].sequence() == cr[10].sequence())

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
    
    // Tests for NoGapAlignmentVector and GappedAlignmentVector
    YMIR_TEST(test_nogap_alignment_vector_no_err())
    YMIR_TEST(test_nogap_alignment_vector_errors())
    YMIR_TEST(test_gapped_alignment_vector())
    YMIR_TEST(test_codon_alignment_vector())

    // Tests for VDJAlignment and VDJAlignmentBuilder
    YMIR_TEST(test_vdj_alignment_nuc_simple_vj())
    YMIR_TEST(test_vdj_alignment_nuc_simple_vdj())
    YMIR_TEST(test_vdj_alignment_nuc_vector_vj())
    YMIR_TEST(test_vdj_alignment_nuc_vector_vdj())

    YMIR_TEST(test_vdj_alignment_aa_simple_vj())
    YMIR_TEST(test_vdj_alignment_aa_simple_vdj())
    YMIR_TEST(test_vdj_alignment_aa_vector_vj())
    YMIR_TEST(test_vdj_alignment_aa_vector_vdj())

    // Tests for clone, clone alignment and clone builder classes.
    YMIR_TEST(test_clone())
    YMIR_TEST(test_clonebuilder_clonealign())
    YMIR_TEST(test_clorep())

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