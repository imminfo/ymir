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


YMIR_TEST_START(test_codon_table)
    YMIR_ASSERT2("G", CodonTable::table().translate("GGC"))
    YMIR_ASSERT2("C", CodonTable::table().translate("TGC"))
    YMIR_ASSERT2("L", CodonTable::table().translate("CTA"))
    YMIR_ASSERT2("T", CodonTable::table().translate("ACG"))
    YMIR_ASSERT2("F", CodonTable::table().translate("TTC"))
    YMIR_ASSERT2("L", CodonTable::table().translate("CTC"))
    YMIR_ASSERT2("N", CodonTable::table().translate("AAC"))
    YMIR_ASSERT2("F", CodonTable::table().translate("TTT"))
    YMIR_ASSERT2("A", CodonTable::table().translate("GCT"))
    YMIR_ASSERT2("R", CodonTable::table().translate("CGG"))
    YMIR_ASSERT2("R", CodonTable::table().translate("CGC"))
    YMIR_ASSERT2("R", CodonTable::table().translate("AGG"))
    YMIR_ASSERT2("L", CodonTable::table().translate("CTT"))
    YMIR_ASSERT2("S", CodonTable::table().translate("TCG"))
    YMIR_ASSERT2("L", CodonTable::table().translate("TTA"))
    YMIR_ASSERT2("V", CodonTable::table().translate("GTT"))
    YMIR_ASSERT2("S", CodonTable::table().translate("AGC"))
    YMIR_ASSERT2("I", CodonTable::table().translate("ATC"))
    YMIR_ASSERT2("A", CodonTable::table().translate("GCA"))
    YMIR_ASSERT2("A", CodonTable::table().translate("GCG"))
    YMIR_ASSERT2("R", CodonTable::table().translate("AGA"))
    YMIR_ASSERT2("V", CodonTable::table().translate("GTG"))
    YMIR_ASSERT2("A", CodonTable::table().translate("GCC"))
    YMIR_ASSERT2("G", CodonTable::table().translate("GGG"))
    YMIR_ASSERT2("T", CodonTable::table().translate("ACA"))
    YMIR_ASSERT2("T", CodonTable::table().translate("ACC"))
    YMIR_ASSERT2("E", CodonTable::table().translate("GAA"))
    YMIR_ASSERT2("K", CodonTable::table().translate("AAA"))
    YMIR_ASSERT2("R", CodonTable::table().translate("CGA"))
    YMIR_ASSERT2("V", CodonTable::table().translate("GTA"))
    YMIR_ASSERT2("P", CodonTable::table().translate("CCT"))
    YMIR_ASSERT2("P", CodonTable::table().translate("CCG"))
    YMIR_ASSERT2("V", CodonTable::table().translate("GTC"))
    YMIR_ASSERT2("R", CodonTable::table().translate("CGT"))
    YMIR_ASSERT2("Y", CodonTable::table().translate("TAT"))
    YMIR_ASSERT2("L", CodonTable::table().translate("CTG"))
    YMIR_ASSERT2("Q", CodonTable::table().translate("CAA"))
    YMIR_ASSERT2("L", CodonTable::table().translate("TTG"))
    YMIR_ASSERT2("H", CodonTable::table().translate("CAC"))
    YMIR_ASSERT2("G", CodonTable::table().translate("GGA"))
    YMIR_ASSERT2("P", CodonTable::table().translate("CCC"))
    YMIR_ASSERT2("N", CodonTable::table().translate("AAT"))
    YMIR_ASSERT2("K", CodonTable::table().translate("AAG"))
    YMIR_ASSERT2("Y", CodonTable::table().translate("TAC"))
    YMIR_ASSERT2("D", CodonTable::table().translate("GAT"))
    YMIR_ASSERT2("Q", CodonTable::table().translate("CAG"))
    YMIR_ASSERT2("S", CodonTable::table().translate("TCC"))
    YMIR_ASSERT2("E", CodonTable::table().translate("GAG"))
    YMIR_ASSERT2("I", CodonTable::table().translate("ATT"))
    YMIR_ASSERT2("I", CodonTable::table().translate("ATA"))
    YMIR_ASSERT2("*", CodonTable::table().translate("TAA"))
    YMIR_ASSERT2("*", CodonTable::table().translate("TGA"))
    YMIR_ASSERT2("H", CodonTable::table().translate("CAT"))
    YMIR_ASSERT2("S", CodonTable::table().translate("TCT"))
    YMIR_ASSERT2("T", CodonTable::table().translate("ACT"))
    YMIR_ASSERT2("*", CodonTable::table().translate("TAG"))
    YMIR_ASSERT2("S", CodonTable::table().translate("TCA"))
    YMIR_ASSERT2("D", CodonTable::table().translate("GAC"))
    YMIR_ASSERT2("M", CodonTable::table().translate("ATG"))
    YMIR_ASSERT2("W", CodonTable::table().translate("TGG"))
    YMIR_ASSERT2("P", CodonTable::table().translate("CCA"))
    YMIR_ASSERT2("G", CodonTable::table().translate("GGT"))
    YMIR_ASSERT2("C", CodonTable::table().translate("TGT"))
    YMIR_ASSERT2("S", CodonTable::table().translate("AGT"))

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

    NaiveCDR3NucleotideAligner nna(genes, VDJAlignerParameters(3));

    YMIR_ASSERT2(nna.alignVar(1, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(1, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(1, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignVar(2, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(2, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(2, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignVar(3, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(3, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(3, "ACGT").len(0), 3)

    YMIR_ASSERT2(nna.alignVar(4, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(4, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(4, "ACGT").len(0), 0)

    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").text_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").len(0), 3)

    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").pattern_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").text_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").pattern_start(0), 3)
    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").len(0), 3)

    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").text_start(0), 4)
    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").pattern_start(0), 5)
    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").len(0), 1)

    YMIR_ASSERT2(nna.alignJoi(4, "ACGG").len(0), 0)

    YMIR_ASSERT2(nna.alignDiv(1, "TTAATAA").size(), 0)

    NaiveCDR3NucleotideAligner nna2(genes, VDJAlignerParameters(2));
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").size(), 2)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").pattern_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").text_start(0), 3)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").len(0), 2)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").pattern_start(1), 1)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").text_start(1), 6)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").len(1), 2)

    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").size(), 3)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(0), 2)

    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(1), 5)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(1), 5)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(1), 2)

    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(2), 5)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(2), 12)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(2), 2)

    YMIR_ASSERT2(nna2.alignDiv(3, "ACTGACGACGGTATCTAC").size(), 5)

YMIR_TEST_END


YMIR_TEST_START(test_cdr3_nuc_aligner)

    // TODO: fix this; tests isn't working due to alignment scores.

    NoGapAlignmentVector vec;

    vector<string> avec1 {"V1", "V2", "V3", "V4"};
    vector<string> svec1 {"ACGTT", "ACGT", "ACG", "TTT"};

    vector<string> avec2 {"J1", "J2", "J3", "J4"};
    vector<string> svec2 {"CGT", "TACGT", "TTCGT", "TTTTT"};

    vector<string> avec3 {"D1", "D2", "D3"};
    vector<string> svec3 {"AA", "AACCTT", "ACT"};
    VDJRecombinationGenes genes("V", avec1, svec1, "J", avec2, svec2, "D", avec3, svec3);

    CDR3NucleotideAligner nna(genes, VDJAlignerParameters(3));

    YMIR_ASSERT2(nna.alignVar(1, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(1, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(1, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignVar(2, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(2, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(2, "ACGT").len(0), 4)

    YMIR_ASSERT2(nna.alignVar(3, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(3, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(3, "ACGT").len(0), 3)

    YMIR_ASSERT2(nna.alignVar(4, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(4, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignVar(4, "ACGT").len(0), 3)
    YMIR_ASSERT(nna.alignVar(4, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(nna.alignVar(4, "ACGT").isMismatch(0, 2))
    YMIR_ASSERT(nna.alignVar(4, "ACGT").isMismatch(0, 3))

    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").text_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").pattern_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(1, "ACGT").len(0), 3)
    YMIR_ASSERT(!nna.alignJoi(1, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(!nna.alignJoi(1, "ACGT").isMismatch(0, 2))
    YMIR_ASSERT(!nna.alignJoi(1, "ACGT").isMismatch(0, 3))

    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").pattern_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(2, "ACGT").len(0), 4)
    YMIR_ASSERT(!nna.alignJoi(2, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(!nna.alignJoi(2, "ACGT").isMismatch(0, 2))
    YMIR_ASSERT(!nna.alignJoi(2, "ACGT").isMismatch(0, 3))
    YMIR_ASSERT(!nna.alignJoi(2, "ACGT").isMismatch(0, 4))

    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").pattern_start(0), 2)
    YMIR_ASSERT2(nna.alignJoi(3, "ACGT").len(0), 4)
    YMIR_ASSERT(nna.alignJoi(3, "ACGT").isMismatch(0, 1))
    YMIR_ASSERT(!nna.alignJoi(3, "ACGT").isMismatch(0, 2))
    YMIR_ASSERT(!nna.alignJoi(3, "ACGT").isMismatch(0, 3))
    YMIR_ASSERT(!nna.alignJoi(3, "ACGT").isMismatch(0, 4))

    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").text_start(0), 1)
    YMIR_ASSERT2(nna.alignJoi(4, "ACGT").pattern_start(0), 2)
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

    YMIR_ASSERT2(nna.alignDiv(1, "TTAATAA").size(), 0)

    CDR3NucleotideAligner nna2(genes, VDJAlignerParameters(2));

    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").size(), 6)

    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").pattern_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").text_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").len(0), 2)
    YMIR_ASSERT(nna2.alignDiv(1, "TTAATAA").isMismatch(0, 1))
    YMIR_ASSERT(nna2.alignDiv(1, "TTAATAA").isMismatch(0, 2))

    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").pattern_start(1), 1)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").text_start(1), 2)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").len(1), 2)
    YMIR_ASSERT(nna2.alignDiv(1, "TTAATAA").isMismatch(1, 1))
    YMIR_ASSERT(!nna2.alignDiv(1, "TTAATAA").isMismatch(1, 2))

    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").pattern_start(5), 1)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").text_start(5), 6)
    YMIR_ASSERT2(nna2.alignDiv(1, "TTAATAA").len(5), 2)
    YMIR_ASSERT(!nna2.alignDiv(1, "TTAATAA").isMismatch(5, 1))
    YMIR_ASSERT(!nna2.alignDiv(1, "TTAATAA").isMismatch(5, 2))

    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").size(), 13 - 2 + 1 + 4)

    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(0), 1)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(0), 6)
    YMIR_ASSERT(!nna2.alignDiv(2, "AAGGTTGGGGGTT").isMismatch(0, 1))
    YMIR_ASSERT(!nna2.alignDiv(2, "AAGGTTGGGGGTT").isMismatch(0, 2))
    YMIR_ASSERT(nna2.alignDiv(2, "AAGGTTGGGGGTT").isMismatch(0, 3))
    YMIR_ASSERT(nna2.alignDiv(2, "AAGGTTGGGGGTT").isMismatch(0, 4))
    YMIR_ASSERT(!nna2.alignDiv(2, "AAGGTTGGGGGTT").isMismatch(0, 5))
    YMIR_ASSERT(!nna2.alignDiv(2, "AAGGTTGGGGGTT").isMismatch(0, 6))

    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").pattern_start(15), 1)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").text_start(15), 12)
    YMIR_ASSERT2(nna2.alignDiv(2, "AAGGTTGGGGGTT").len(15), 2)
    YMIR_ASSERT(nna2.alignDiv(2, "AAGGTTGGGGGTT").isMismatch(15, 1))
    YMIR_ASSERT(nna2.alignDiv(2, "AAGGTTGGGGGTT").isMismatch(15, 2))

    YMIR_ASSERT2(nna2.alignDiv(3, "ACTGACGACGGTATCTAC").size(), 18 - 2 + 1 + 1)

    YMIR_ASSERT2(nna2.alignDiv(3, "ACTGACGACGGTATCTAC").pattern_start(17), 1)
    YMIR_ASSERT2(nna2.alignDiv(3, "ACTGACGACGGTATCTAC").text_start(17), 17)
    YMIR_ASSERT2(nna2.alignDiv(3, "ACTGACGACGGTATCTAC").len(17), 2)

YMIR_TEST_END


YMIR_TEST_START(test_cdr3_aa_aligner)

    CodonAlignmentVector vec;

    vector<string> avec1 {"V1", "V2", "V3", "V4", "V5"};
    vector<string> svec1 {"TCGTT", "TCGT", "TCGTTC", "TGGG", "TCGG"};

    vector<string> avec2 {"J1", "J2", "J3", "J4", "J5", "J6"};
    vector<string> svec2 {"CTTTA", "CCTT", "AGCCTG", "AGGCTG", "AATT", "TTT"};

    vector<string> avec3 {"D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10"};
    vector<string> svec3 {"AAA", "AACCACT", "ACT", "TGC", "GGCA", "GGTT", "GGGG", "CCCC", "CCC", "GGG"};
    VDJRecombinationGenes genes("V", avec1, svec1, "J", avec2, svec2, "D", avec3, svec3);

    CDR3AminoAcidAligner aligner(genes, VDJAlignerParameters(3,
                                                             VDJAlignmentEventScore(AlignmentEventScore(1, -1, 1),
                                                                                    AlignmentEventScore(1, -1, 1),
                                                                                    AlignmentEventScore(1, -1, 1)),
                                                             VDJAlignmentScoreThreshold(2, 2, 2)));

    // {'S', "TCT"}, {'S', "TCC"}, {'S', "TCA"}, {'S', "TCG"}, {'S', "AGT"}, {'S', "AGC"},
    // {'L', "TTA"}, {'L', "TTG"}, {'L', "CTT"}, {'L', "CTC"}, {'L', "CTA"}, {'L', "CTG"},
    // V1: TCGTT
    auto alignment = aligner.alignVar(1, "SL");
    YMIR_ASSERT2(alignment.id(0), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 5)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, false, false, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 4), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 5), compute_codon_hash({true, true, false, false, false, false}, 0))

    alignment = aligner.alignVar(2, "SL");
    YMIR_ASSERT2(alignment.id(0), 2)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 4)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, false, false, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 4), compute_codon_hash({true, true, false, false, false, false}, 0))

    alignment = aligner.alignVar(3, "SL");
    YMIR_ASSERT2(alignment.id(0), 3)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 5)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, false, false, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 4), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 5), compute_codon_hash({true, true, false, false, false, false}, 0))

    // {'C', "TGT"}, {'C', "TGC"},
    // TG.GG
    // TG.T
    alignment = aligner.alignVar(4, "CL");
    YMIR_ASSERT2(alignment.id(0), 4)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 2)

    // TCG.C
    // TCG.G
    alignment = aligner.alignVar(5, "SL");
    YMIR_ASSERT2(alignment.id(0), 5)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 3)

    // {'S', "TCT"}, {'S', "TCC"}, {'S', "TCA"}, {'S', "TCG"}, {'S', "AGT"}, {'S', "AGC"},
    // {'L', "TTA"}, {'L', "TTG"}, {'L', "CTT"}, {'L', "CTC"}, {'L', "CTA"}, {'L', "CTG"},
    // J1: CT.TTA
    //    TCT.TTA
    alignment = aligner.alignJoi(1, "SL");
    YMIR_ASSERT2(alignment.id(0), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 2)
    YMIR_ASSERT2(alignment.len(0), 5)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, false, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, false, false, false, true, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({true, false, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 4), compute_codon_hash({true, false, false, false, true, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 5), compute_codon_hash({true, false, false, false, true, false}, 0))

    // J2: C.CTT
    //   TCC.CTT
    //   AGC.TTA
    alignment = aligner.alignJoi(2, "SL");
    YMIR_ASSERT2(alignment.id(0), 2)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 3)
    YMIR_ASSERT2(alignment.len(0), 4)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({false, true, false, false, false, true}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({false, false, true, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, false, true, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 4), compute_codon_hash({false, false, true, false, false, false}, 0))

    // J3: AGC.CTG
    //     AGC.CTG
    alignment = aligner.alignJoi(3, "SL");
    YMIR_ASSERT2(alignment.id(0), 3)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 6)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({false, false, false, false, false, true}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({false, false, false, false, false, true}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, true, false, false, false, true}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 4), compute_codon_hash({false, false, false, false, false, true}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 6), compute_codon_hash({false, true, false, false, false, true}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 5), compute_codon_hash({false, true, false, false, false, true}, 0))

    // J4: AGG.CTG
    //     TCG.CTG
    alignment = aligner.alignJoi(4, "SL");
    YMIR_ASSERT2(alignment.id(0), 4)
    YMIR_ASSERT2(alignment.pattern_start(0), 3)
    YMIR_ASSERT2(alignment.text_start(0), 3)
    YMIR_ASSERT2(alignment.len(0), 4)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({false, false, false, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({false, false, false, false, false, true}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, true, false, false, false, true}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 4), compute_codon_hash({false, true, false, false, false, true}, 0))

    // J5: AATT
    // S: TCT TCC TCA TCG AGT AGC
    // F: TTT TTC
    alignment = aligner.alignJoi(5, "SF");
    YMIR_ASSERT2(alignment.size(), 1)
    YMIR_ASSERT2(alignment.id(0), 5)
    YMIR_ASSERT2(alignment.pattern_start(0), 3)
    YMIR_ASSERT2(alignment.text_start(0), 5)
    YMIR_ASSERT2(alignment.len(0), 2)

    alignment = aligner.alignJoi(6, "ASF");
    YMIR_ASSERT2(alignment.size(), 1)
    YMIR_ASSERT2(alignment.id(0), 6)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 7)
    YMIR_ASSERT2(alignment.len(0), 3)

    // "D1", "D2", "D3"
    // "AAA", "AACCACT", "ACT", "TGC"
    // one aminoacid
    // {'K', "AAA"}, {'K', "AAG"}
    // {'Q', "CAA"}, {'Q', "CAG"}
    // {'H', "CAT"}, {'H', "CAC"}
    // {'T', "ACT"}, {'T', "ACC"}, {'T', "ACA"}, {'T', "ACG"}
    // {'S', "TCT"}, {'S', "TCC"}, {'S', "TCA"}, {'S', "TCG"}, {'S', "AGT"}, {'S', "AGC"}
    alignment = aligner.alignDiv(1, "Q");
    YMIR_ASSERT2(alignment.size(), 0)

    alignment = aligner.alignDiv(1, "K");
    YMIR_ASSERT2(alignment.size(), 1)
    YMIR_ASSERT2(alignment.id(0), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 3)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({true, false, false, false, false, false}, 0))

    alignment = aligner.alignDiv(1, "HK");
    YMIR_ASSERT2(alignment.size(), 1)
    YMIR_ASSERT2(alignment.id(0), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 4)
    YMIR_ASSERT2(alignment.len(0), 3)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({true, false, false, false, false, false}, 0))


    alignment = aligner.alignDiv(2, "T");
    YMIR_ASSERT2(alignment.size(), 2)

    YMIR_ASSERT2(alignment.id(0), 2)
    YMIR_ASSERT2(alignment.pattern_start(0), 2)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 3)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, true, false, false, false, false}, 0))

    YMIR_ASSERT2(alignment.id(1), 2)
    YMIR_ASSERT2(alignment.pattern_start(1), 5)
    YMIR_ASSERT2(alignment.text_start(1), 1)
    YMIR_ASSERT2(alignment.len(1), 3)
    YMIR_ASSERTCOD(alignment.getCodon(1, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(1, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(1, 3), compute_codon_hash({true, false, false, false, false, false}, 0))


    alignment = aligner.alignDiv(3, "T");
    YMIR_ASSERT2(alignment.size(), 1)
    YMIR_ASSERT2(alignment.id(0), 3)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 3)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({true, false, false, false, false, false}, 0))


    alignment = aligner.alignDiv(4, "S");
    YMIR_ASSERT2(alignment.size(), 0)

    // neighbour aminoacids
    // "D1", "D2", "D3"
    // "AAA", "AACCACT", "ACT"
    //
    // {'Q', "CAA"}, {'Q', "CAG"}
    // {'R', "CGT"}, {'R', "CGC"}, {'R', "CGA"}, {'R', "CGG"}, {'R', "AGA"}, {'R', "AGG"}
    // {'K', "AAA"}, {'K', "AAG"}
    alignment = aligner.alignDiv(1, "QRK");
    YMIR_ASSERT2(alignment.size(), 3)

    YMIR_ASSERT2(alignment.id(0), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 2)
    YMIR_ASSERT2(alignment.len(0), 3)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, false, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, false, false, false, true, true}, 0))

    YMIR_ASSERT2(alignment.id(1), 1)
    YMIR_ASSERT2(alignment.pattern_start(1), 1)
    YMIR_ASSERT2(alignment.text_start(1), 6)
    YMIR_ASSERT2(alignment.len(1), 3)
    YMIR_ASSERTCOD(alignment.getCodon(1, 1), compute_codon_hash({false, false, true, false, true, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(1, 2), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(1, 3), compute_codon_hash({true, true, false, false, false, false}, 0))

    YMIR_ASSERT2(alignment.id(2), 1)
    YMIR_ASSERT2(alignment.pattern_start(2), 1)
    YMIR_ASSERT2(alignment.text_start(2), 7)
    YMIR_ASSERT2(alignment.len(2), 3)
    YMIR_ASSERTCOD(alignment.getCodon(2, 1), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(2, 2), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(2, 3), compute_codon_hash({true, false, false, false, false, false}, 0))

    // three aminoacids
    // "D1", "D2", "D3"
    // "AAA", "AACCACT", "ACT"
    //
    // {'N', "AAT"}, {'N', "AAC"},
    // {'H', "CAT"}, {'H', "CAC"},
    // {'F', "TTT"}, {'F', "TTC"},
    // {'N', "AAT"}, {'N', "AAC"},
    // {'T', "ACT"}, {'T', "ACC"}, {'T', "ACA"}, {'T', "ACG"}
    alignment = aligner.alignDiv(2, "NHFNT");
    YMIR_ASSERT2(alignment.size(), 4)

    YMIR_ASSERT2(alignment.id(0), 2)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 7)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 4), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 5), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 6), compute_codon_hash({false, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 7), compute_codon_hash({true, true, false, false, false, false}, 0))

    YMIR_ASSERT2(alignment.id(1), 2)
    YMIR_ASSERT2(alignment.pattern_start(1), 4)
    YMIR_ASSERT2(alignment.text_start(1), 12)
    YMIR_ASSERT2(alignment.len(1), 4)
    YMIR_ASSERTCOD(alignment.getCodon(1, 1), compute_codon_hash({false, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(1, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(1, 3), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(1, 4), compute_codon_hash({true, false, false, false, false, false}, 0))

    YMIR_ASSERT2(alignment.id(2), 2)
    YMIR_ASSERT2(alignment.pattern_start(2), 1)
    YMIR_ASSERT2(alignment.text_start(2), 10)
    YMIR_ASSERT2(alignment.len(2), 3)
    YMIR_ASSERTCOD(alignment.getCodon(2, 1), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(2, 2), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(2, 3), compute_codon_hash({false, true, false, false, false, false}, 0))

    YMIR_ASSERT2(alignment.id(3), 2)
    YMIR_ASSERT2(alignment.pattern_start(3), 2)
    YMIR_ASSERT2(alignment.text_start(3), 13)
    YMIR_ASSERT2(alignment.len(3), 3)
    YMIR_ASSERTCOD(alignment.getCodon(3, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(3, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(3, 3), compute_codon_hash({false, true, false, false, false, false}, 0))


    alignment = aligner.alignDiv(5, "CAGVQF");
    YMIR_ASSERT2(alignment.size(), 3)

    YMIR_ASSERT2(alignment.id(0), 5)
    YMIR_ASSERT2(alignment.pattern_start(0), 2)
    YMIR_ASSERT2(alignment.text_start(0), 4)
    YMIR_ASSERT2(alignment.len(0), 3)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, false, true, false, false, false}, 0))

    YMIR_ASSERT2(alignment.id(1), 5)
    YMIR_ASSERT2(alignment.pattern_start(1), 1)
    YMIR_ASSERT2(alignment.text_start(1), 7)
    YMIR_ASSERT2(alignment.len(1), 3)
    YMIR_ASSERTCOD(alignment.getCodon(1, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(1, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(1, 3), compute_codon_hash({false, true, false, false, false, false}, 0))

    YMIR_ASSERT2(alignment.id(2), 5)
    YMIR_ASSERT2(alignment.pattern_start(2), 2)
    YMIR_ASSERT2(alignment.text_start(2), 12)
    YMIR_ASSERT2(alignment.len(2), 3)
    YMIR_ASSERTCOD(alignment.getCodon(2, 1), compute_codon_hash({false, false, false, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(2, 2), compute_codon_hash({true, true, false, false, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(2, 3), compute_codon_hash({true, true, false, false, false, false}, 0))


    alignment = aligner.alignDiv(7, "CAGF");
    YMIR_ASSERT2(alignment.size(), 3)


    alignment = aligner.alignDiv(8, "CALF");
    YMIR_ASSERT2(alignment.size(), 2)

    alignment = aligner.alignDiv(9, "KPF");
    YMIR_ASSERT2(alignment.size(), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 4)
    YMIR_ASSERT2(alignment.len(0), 3)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, true, false, false, false, false}, 0))

    alignment = aligner.alignDiv(10, "NGW");
    YMIR_ASSERT2(alignment.size(), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 4)
    YMIR_ASSERT2(alignment.len(0), 3)
    YMIR_ASSERTCOD(alignment.getCodon(0, 1), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 2), compute_codon_hash({true, true, true, true, false, false}, 0))
    YMIR_ASSERTCOD(alignment.getCodon(0, 3), compute_codon_hash({false, false, false, true, false, false}, 0))

YMIR_TEST_END


YMIR_TEST_START(test_sw_alignment_matrix)

    /*
          1  2  3  4  5  6
       0  0  0  0  0  0  0  0
     1 0 *1  2  0
     2 0  0  3  0
     3 0  0  0 *1--2\ 2
     4 0  0  0  0  0 |3  3
     5 0  0  0  0  0 |4--5

     */
    sequence_t alpha, beta;
    alpha = "123456";
    beta = "12345";

    SWAlignmentMatrix mat(10, beta, alpha);

    mat.setStart(1, 1);
    mat.setStart(3, 3);

    mat.score(1, 1) = 1;
    mat.score(1, 2) = 2;
    mat.score(2, 2) = 3;
    mat.score(3, 3) = 1;
    mat.score(3, 4) = 2;
    mat.score(3, 5) = 2;
    mat.score(4, 5) = 3;
    mat.score(5, 5) = 4;
    mat.score(5, 6) = 5;
    mat.score(4, 6) = 3;

    GappedAlignmentVector vec;
    YMIR_ASSERT2(mat.getBestAlignment(&vec, beta, alpha), 5);
    YMIR_ASSERT2(vec.id(0), 10)
    YMIR_ASSERT2(vec.pattern_start(0), 3)
    YMIR_ASSERT2(vec.text_start(0), 3)
    YMIR_ASSERT2(vec.len(0), 5)

    YMIR_ASSERT(vec.isMatch(0, 1))
    YMIR_ASSERT(vec.isIns(0, 2))
    YMIR_ASSERT(vec.isMatch(0, 3))
    YMIR_ASSERT(vec.isDel(0, 4))
    YMIR_ASSERT(vec.isIns(0, 5))

YMIR_TEST_END


YMIR_TEST_START(test_sw_aligner)

    YMIR_ASSERT(false)

    vector<string> avec1 {"V1", "V2", "V3", "V4"};
    vector<string> svec1 {"ACGTTGGC", "ACGTTAGCTA", "ACG", "TTT"};

    vector<string> avec2 {"J1", "J2", "J3", "J4"};
    vector<string> svec2 {"CGT", "TACGT", "TTCGT", "TTTTT"};

    vector<string> avec3 {"D1", "D2", "D3"};
    vector<string> svec3 {"AA", "AACCTT", "ACT"};

    VDJRecombinationGenes genes("V", avec1, svec1, "J", avec2, svec2, "D", avec3, svec3);
//
//    SmithWatermanNoGapAligner swa(genes, VDJAlignerParameters(1,
//                                                              3,
//                                                              AlignmentEventScore(1, -1, -3),
//                                                              AlignmentEventScore(2, -1, -1),
//                                                              AlignmentEventScore(3, -2, -3)));


YMIR_TEST_END


YMIR_TEST_START(test_swng_alignment_matrix)

    /*
     0  0  0  0  0  0  0  0
     0 *1  0  0
     0  0  2 *1
     0  0  0  3  2
     0 *1  0  0  0  3     0
     0  0  2  0  0  0  4  3

     */
    sequence_t alpha, beta;
    alpha = "1123657";
    beta =   "12345";

    SWNGAlignmentMatrix mat(10, beta, alpha);

    mat.score(1, 1) = 1;
    mat.score(2, 2) = 2;
    mat.score(3, 3) = 3;
    mat.score(2, 3) = 1;
    mat.score(3, 4) = 2;
    mat.score(4, 5) = 3;
    mat.score(5, 6) = 4;
    mat.score(4, 1) = 1;
    mat.score(5, 2) = 2;
    mat.score(5, 7) = 3;

    mat.setStart(1, 1);
    mat.setStart(2, 3);
    mat.setStart(4, 1);

    NoGapAlignmentVector vec;
    YMIR_ASSERT2(mat.getBestAlignment(&vec, beta, alpha), 4);
    YMIR_ASSERT2(vec.id(0), 10)
    YMIR_ASSERT2(vec.pattern_start(0), 1)
    YMIR_ASSERT2(vec.text_start(0), 2)
    YMIR_ASSERT2(vec.len(0), 5)

    YMIR_ASSERT(!vec.isMismatch(0, 1))
    YMIR_ASSERT(!vec.isMismatch(0, 2))
    YMIR_ASSERT(!vec.isMismatch(0, 3))
    YMIR_ASSERT(vec.isMismatch(0, 4))
    YMIR_ASSERT(!vec.isMismatch(0, 5))

YMIR_TEST_END


YMIR_TEST_START(test_swng_aligner)

    vector<string> avec1 {"V1", "V2", "V3", "V4"};
    vector<string> svec1 {"ACGTTGGGATC", "ACGT", "ACG", "TTT"};

    vector<string> avec2 {"J1", "J2", "J3"};
    vector<string> svec2 {"TTT", "ACGTTGGGATC", "ACGT"};

    vector<string> avec3 {"D1", "D2", "D3"};
    vector<string> svec3 {"AA", "AACCTT", "ACT"};

    VDJRecombinationGenes genes("V", avec1, svec1, "J", avec2, svec2, "D", avec3, svec3);

    SmithWatermanNoGapAligner swnga(genes, VDJAlignerParameters(2,
                                                                VDJAlignmentEventScore(AlignmentEventScore(3, -1, -3),
                                                                                       AlignmentEventScore(2, -1, -1),
                                                                                       AlignmentEventScore(3, -2, -3))));

    sequence_t pattern;
    NoGapAlignmentVector alignment;

    // V: ACGTTGGGATC
    // P:    ATGGGTT
    pattern = "ATGGGTT";
    alignment = swnga.alignVar(1, pattern);
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 4)
    YMIR_ASSERT2(alignment.len(0), 7)
    YMIR_ASSERT2(alignment.id(0), 1)
    YMIR_ASSERT(alignment.isMismatch(0, 1))
    YMIR_ASSERT(!alignment.isMismatch(0, 2))
    YMIR_ASSERT(!alignment.isMismatch(0, 3))
    YMIR_ASSERT(!alignment.isMismatch(0, 4))
    YMIR_ASSERT(!alignment.isMismatch(0, 5))
    YMIR_ASSERT(alignment.isMismatch(0, 6))
    YMIR_ASSERT(!alignment.isMismatch(0, 7))

    alignment = swnga.alignJoi(2, pattern);
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 4)
    YMIR_ASSERT2(alignment.len(0), 7)
    YMIR_ASSERT2(alignment.id(0), 2)
    YMIR_ASSERT(alignment.isMismatch(0, 1))
    YMIR_ASSERT(!alignment.isMismatch(0, 2))
    YMIR_ASSERT(!alignment.isMismatch(0, 3))
    YMIR_ASSERT(!alignment.isMismatch(0, 4))
    YMIR_ASSERT(!alignment.isMismatch(0, 5))
    YMIR_ASSERT(alignment.isMismatch(0, 6))
    YMIR_ASSERT(!alignment.isMismatch(0, 7))

    // V: ACGTTGGGATC
    // P:       TGATCTT
    pattern = "TGATCTT";
    alignment = swnga.alignVar(1, pattern);
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 7)
    YMIR_ASSERT2(alignment.len(0), 5)
    YMIR_ASSERT2(alignment.id(0), 1)
    YMIR_ASSERT(alignment.isMismatch(0, 1))
    YMIR_ASSERT(!alignment.isMismatch(0, 2))
    YMIR_ASSERT(!alignment.isMismatch(0, 3))
    YMIR_ASSERT(!alignment.isMismatch(0, 4))
    YMIR_ASSERT(!alignment.isMismatch(0, 5))

    alignment = swnga.alignJoi(2, pattern);
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.pattern_start(0), 7)
    YMIR_ASSERT2(alignment.len(0), 5)
    YMIR_ASSERT2(alignment.id(0), 2)
    YMIR_ASSERT(alignment.isMismatch(0, 1))
    YMIR_ASSERT(!alignment.isMismatch(0, 2))
    YMIR_ASSERT(!alignment.isMismatch(0, 3))
    YMIR_ASSERT(!alignment.isMismatch(0, 4))
    YMIR_ASSERT(!alignment.isMismatch(0, 5))

    // V:    ACGT
    // P: GGAACTTAGGG
    pattern = "GGAACTTAGGG";
    alignment = swnga.alignVar(2, pattern);
    YMIR_ASSERT2(alignment.text_start(0), 4)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 4)
    YMIR_ASSERT2(alignment.id(0), 2)
    YMIR_ASSERT(!alignment.isMismatch(0, 1))
    YMIR_ASSERT(!alignment.isMismatch(0, 2))
    YMIR_ASSERT(alignment.isMismatch(0, 3))
    YMIR_ASSERT(!alignment.isMismatch(0, 4))

    alignment = swnga.alignJoi(3, pattern);
    YMIR_ASSERT2(alignment.text_start(0), 4)
    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 4)
    YMIR_ASSERT2(alignment.id(0), 3)
    YMIR_ASSERT(!alignment.isMismatch(0, 1))
    YMIR_ASSERT(!alignment.isMismatch(0, 2))
    YMIR_ASSERT(alignment.isMismatch(0, 3))
    YMIR_ASSERT(!alignment.isMismatch(0, 4))

    alignment = swnga.alignDiv(1, "TTAATAA");
    YMIR_ASSERT2(alignment.size(), 6)

    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 2)
    YMIR_ASSERT(alignment.isMismatch(0, 1))
    YMIR_ASSERT(alignment.isMismatch(0, 2))

    YMIR_ASSERT2(alignment.pattern_start(1), 1)
    YMIR_ASSERT2(alignment.text_start(1), 2)
    YMIR_ASSERT2(alignment.len(1), 2)
    YMIR_ASSERT(alignment.isMismatch(1, 1))
    YMIR_ASSERT(!alignment.isMismatch(1, 2))

    YMIR_ASSERT2(alignment.pattern_start(5), 1)
    YMIR_ASSERT2(alignment.text_start(5), 6)
    YMIR_ASSERT2(alignment.len(5), 2)
    YMIR_ASSERT(!alignment.isMismatch(5, 1))
    YMIR_ASSERT(!alignment.isMismatch(5, 2))

    alignment = swnga.alignDiv(2, "AAGGTTGGGGGTT");
    YMIR_ASSERT2(alignment.size(), 13 - 2 + 1 + 4)

    YMIR_ASSERT2(alignment.pattern_start(0), 1)
    YMIR_ASSERT2(alignment.text_start(0), 1)
    YMIR_ASSERT2(alignment.len(0), 6)
    YMIR_ASSERT(!alignment.isMismatch(0, 1))
    YMIR_ASSERT(!alignment.isMismatch(0, 2))
    YMIR_ASSERT(alignment.isMismatch(0, 3))
    YMIR_ASSERT(alignment.isMismatch(0, 4))
    YMIR_ASSERT(!alignment.isMismatch(0, 5))
    YMIR_ASSERT(!alignment.isMismatch(0, 6))

    YMIR_ASSERT2(alignment.pattern_start(15), 1)
    YMIR_ASSERT2(alignment.text_start(15), 12)
    YMIR_ASSERT2(alignment.len(15), 2)
    YMIR_ASSERT(alignment.isMismatch(15, 1))
    YMIR_ASSERT(alignment.isMismatch(15, 2))

    alignment = swnga.alignDiv(3, "ACTGACGACGGTATCTAC");
    YMIR_ASSERT2(alignment.size(), 18 - 2 + 1 + 1)

    YMIR_ASSERT2(alignment.pattern_start(17), 1)
    YMIR_ASSERT2(alignment.text_start(17), 17)
    YMIR_ASSERT2(alignment.len(17), 2)

YMIR_TEST_END


YMIR_TEST_START(test_errcorr_aligner)

    YMIR_ASSERT(false)

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
    // Test for CodonTable.
//    YMIR_TEST(test_codon_table())

    // Tests for sequences aligners.
//    YMIR_TEST(test_naive_cdr3_nuc_aligner())
//    YMIR_TEST(test_cdr3_nuc_aligner())
    YMIR_TEST(test_cdr3_aa_aligner())
//    YMIR_TEST(test_sw_alignment_matrix())
//    YMIR_TEST(test_sw_aligner())
//    YMIR_TEST(test_swng_alignment_matrix())
//    YMIR_TEST(test_swng_aligner())

    // Error corrector test
//    YMIR_TEST(test_errcorr_aligner())

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