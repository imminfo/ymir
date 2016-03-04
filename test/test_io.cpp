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


YMIR_TEST_START(test_writer)

     RepertoireWriter writer;

     vector<string> alvec1;
     vector<string> seqvec1;
     alvec1.push_back("Vseg1");
     alvec1.push_back("Vseg2");
     alvec1.push_back("Vseg3");
     seqvec1.push_back("CCCG");
     seqvec1.push_back("GGG");
     seqvec1.push_back("AGGCGAG");

     vector<string> alvec2;
     vector<string> seqvec2;
     alvec2.push_back("Jseg1");
     alvec2.push_back("Jseg2");
     alvec2.push_back("Jseg3");
     seqvec2.push_back("CCGTTT");
     seqvec2.push_back("ATTTGG");
     seqvec2.push_back("AGGTTT");

     VDJRecombinationGenes genes("VA", alvec1, seqvec1, "JA", alvec2, seqvec2);

     ClonotypeBuilder cl_builder;
     // CCCG.AC.GGTTT
     cl_builder.setSequence("CCCGACGGTTT")
             .setNucleotideSeq()
             .setRecombination(VJ_RECOMB)
             .addVarAlignment(1, 1, 1, 4)
             .addVarAlignment(3, 4, 3, 4)
             .addJoiAlignment(1, 2, 6, 5)
             .addJoiAlignment(2, 2, 9, 3)
             .addJoiAlignment(3, 2, 7, 5);
     Clonotype clonotype = cl_builder.buildClonotype();
     vector<Clonotype> vec;
     vec.push_back(clonotype);
     Cloneset cloneset(vec);

     YMIR_ASSERT(writer.write(TEST_DATA_FOLDER + "../out.txt", cloneset, genes))

//     TODO: write a parser to test writer's output

YMIR_TEST_END


YMIR_TEST_START(test_parser_vj)

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

     YMIR_ASSERT(cr.size() == 30)
     YMIR_ASSERT(cr[0].sequence() == "TGTGCAGCAAGTACCCCCTTAAGCTGGTGGTACTAGCTATGGAAAGCTGACATTT")
     YMIR_ASSERT(cr[0].recombination() == VJ_RECOMB)
     YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(0)].allele == "TRAV13-1")
     YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(1)].allele == "TRAV13-2")
     YMIR_ASSERT(cr[0].nVar() == 2)
     YMIR_ASSERT(vdj_genes.J()[cr[0].getJoi(0)].allele == "TRAJ52")
     YMIR_ASSERT(cr[0].nJoi() == 1)

     YMIR_ASSERT(cr[2].sequence() == "TGTGCAACTCTTAGCAGGGATGAACACAGGCTTTCAGAAACTTGTATTT")
     YMIR_ASSERT(cr[2].recombination() != VDJ_RECOMB)
     YMIR_ASSERT(vdj_genes.V()[cr[2].getVar(0)].allele == "TRAV12-3")
     YMIR_ASSERT(cr[2].nVar() == 1)
     YMIR_ASSERT(vdj_genes.J()[cr[2].getJoi(0)].allele == "TRAJ8")
     YMIR_ASSERT(vdj_genes.J()[cr[2].getJoi(1)].allele == "TRAJ18")
     YMIR_ASSERT(cr[2].nJoi() == 2)

YMIR_TEST_END


YMIR_TEST_START(test_parser_vj_by_block)

    NaiveNucParser parser;

    bool V_err, J_err;
    VDJRecombinationGenes vdj_genes("Vgene", TEST_DATA_FOLDER + "vgene.real.txt"
            , "Jgene", TEST_DATA_FOLDER + "jgene.real.txt", &V_err, &J_err);
    YMIR_ASSERT(V_err)
    YMIR_ASSERT(J_err)

    Cloneset cr;
    YMIR_ASSERT(parser.open(TEST_DATA_FOLDER + "ymir.alpha.txt",
                                    vdj_genes,
                                    NUCLEOTIDE,
                                    VJ_RECOMB,
                                    AlignmentColumnOptions()
                                            .setV(AlignmentColumnOptions::USE_PROVIDED)
                                            .setJ(AlignmentColumnOptions::USE_PROVIDED)))

    YMIR_ASSERT(parser.parse(&cr, 1))
    YMIR_ASSERT2(cr.size(), 1)
    YMIR_ASSERT(cr[0].sequence() == "TGTGCAGCAAGTACCCCCTTAAGCTGGTGGTACTAGCTATGGAAAGCTGACATTT")
    YMIR_ASSERT(cr[0].recombination() == VJ_RECOMB)
    YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(0)].allele == "TRAV13-1")
    YMIR_ASSERT(vdj_genes.V()[cr[0].getVar(1)].allele == "TRAV13-2")
    YMIR_ASSERT(cr[0].nVar() == 2)
    YMIR_ASSERT(vdj_genes.J()[cr[0].getJoi(0)].allele == "TRAJ52")
    YMIR_ASSERT(cr[0].nJoi() == 1)

    YMIR_ASSERT(parser.parse(&cr, 2))
    YMIR_ASSERT2(cr.size(), 2)
    YMIR_ASSERT2(cr[1].sequence(),"TGTGCAACTCTTAGCAGGGATGAACACAGGCTTTCAGAAACTTGTATTT")
    YMIR_ASSERT(cr[1].recombination() != VDJ_RECOMB)
    YMIR_ASSERT2(cr[1].nVar(), 1)
    YMIR_ASSERT2(cr[1].nJoi(), 2)
    YMIR_ASSERT(vdj_genes.V()[cr[1].getVar(0)].allele == "TRAV12-3")
    YMIR_ASSERT(vdj_genes.J()[cr[1].getJoi(0)].allele == "TRAJ8")
    YMIR_ASSERT(vdj_genes.J()[cr[1].getJoi(1)].allele == "TRAJ18")

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

    NaiveNucParser parser;

    Cloneset cr;
    YMIR_ASSERT(parser.openAndParse(TEST_DATA_FOLDER + "ymir.beta.txt",
                                     &cr,
                                     genes,
                                     NUCLEOTIDE,
                                     VDJ_RECOMB,
                                     AlignmentColumnOptions()
                                      .setV(AlignmentColumnOptions::USE_PROVIDED)
                                      .setJ(AlignmentColumnOptions::USE_PROVIDED)
                                      .setD(AlignmentColumnOptions::OVERWRITE),
                                    VDJAlignerParameters(1, 3)))

    YMIR_ASSERT2(cr.size(), 1)
    YMIR_ASSERT2(cr[0].sequence(), "CCCGACGGTTT")
    YMIR_ASSERT2(genes.V()[cr[0].getJoi(0)].allele, "Vseg1")
    YMIR_ASSERT2(genes.J()[cr[0].getJoi(0)].allele, "Jseg1")

    // CCCGACGGTTT
    // .......GTTT
    // ...ACCGGT
    // ACCGGT
    //CCCGGAC
    // CCCGGAC
    // ...CCCGGAC
    YMIR_ASSERT2( (int) cr[0].nDiv(), 3)
    YMIR_ASSERT2(cr[0].getDiv(0), 1);
    YMIR_ASSERT2(cr[0].getDiv(1), 2);
    YMIR_ASSERT2(cr[0].getDiv(2), 3);
    YMIR_ASSERT2( (int) cr[0].numDivAlignments(0), 1)
    YMIR_ASSERT2( (int) cr[0].numDivAlignments(1), 2)
    YMIR_ASSERT2( (int) cr[0].numDivAlignments(2), 3)
    YMIR_ASSERT2(cr[0].getDivSeqStart(0, 0), 8)
    YMIR_ASSERT2(cr[0].getDivGeneStart(0, 0), 1)
    YMIR_ASSERT2(cr[0].getDivLen(0, 0), 4)
    YMIR_ASSERT2(cr[0].getDivSeqStart(1, 0), 2)
    YMIR_ASSERT2(cr[0].getDivGeneStart(1, 0), 2)
    YMIR_ASSERT2(cr[0].getDivLen(1, 0), 3)
    YMIR_ASSERT2(cr[0].getDivSeqStart(1, 1), 6)
    YMIR_ASSERT2(cr[0].getDivGeneStart(1, 1), 3)
    YMIR_ASSERT2(cr[0].getDivLen(1, 1), 4)

    YMIR_ASSERT(parser.openAndParse(TEST_DATA_FOLDER + "ymir.beta.txt",
                                    &cr,
                                    genes,
                                    NUCLEOTIDE,
                                    VDJ_RECOMB,
                                    AlignmentColumnOptions()
                                            .setV(AlignmentColumnOptions::USE_PROVIDED)
                                            .setJ(AlignmentColumnOptions::USE_PROVIDED)
                                            .setD(AlignmentColumnOptions::OVERWRITE),
                                    VDJAlignerParameters(1, 2)))
//     this is if default gene len is equal to 2
    YMIR_ASSERT2( (int) cr[0].numDivAlignments(0), 3)
    YMIR_ASSERT2( (int) cr[0].numDivAlignments(1), 4)
    YMIR_ASSERT2( (int) cr[0].numDivAlignments(2), 5)

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

    // Test for the repertoire parser and writer
    YMIR_TEST(test_parser_vj())
    YMIR_TEST(test_parser_vj_by_block())
    YMIR_TEST(test_parser_vdj_with_d_alignment())
    YMIR_TEST(test_writer())

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