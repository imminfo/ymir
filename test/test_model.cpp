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

    YMIR_ASSERT(false)

YMIR_TEST_END


YMIR_TEST_START(test_model_vj_maag)

    // TGTGCTCTTGGGGAACTTTCGGAGTGGCTCTAGCAACACAGGCAAACTAATCTTT	CALGELSEW~SSNTGKLIF	TRAV13-2		TRAJ37	1|1|5		0|25|31

    ModelParameterVector mvec = make_test_events_vj();

    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vj_model/");
    YMIR_ASSERT(model.status())

     ParserNuc parser(new NaiveCDR3NucleotideAligner());

     bool V_err, J_err;
     VDJRecombinationGenes vj_genes("Vgene", TEST_DATA_FOLDER + "vgene.real.txt",
                                     "Jgene", TEST_DATA_FOLDER + "jgene.real.txt",
                                     &V_err, &J_err);
     YMIR_ASSERT(V_err)
     YMIR_ASSERT(J_err)

     ClonesetNuc cloneset;
    YMIR_ASSERT(parser.openAndParse(TEST_DATA_FOLDER + "ymir.alpha2.txt",
                                    &cloneset,
                                    vj_genes,
                                    VJ_RECOMB,
                                    AlignmentColumnOptions()
                                            .setV(AlignmentColumnOptions::USE_PROVIDED)
                                            .setJ(AlignmentColumnOptions::USE_PROVIDED)))

     MAAGnuc maag = model.buildGraphs(cloneset, SAVE_METADATA, NO_ERRORS)[1];

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


YMIR_TEST_START(test_model_vdj_maag)

    ProbabilisticAssemblingModel model(TEST_DATA_FOLDER + "test_vdj_model/");
    YMIR_ASSERT(model.status())

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
    
    // Tests for probabilistic assembling model (PAM) reading / writing files.
    YMIR_TEST(test_model_vj_file())
    YMIR_TEST(test_model_vdj_file())
    YMIR_TEST(test_model_vj_save_load())
    YMIR_TEST(test_model_vdj_save_load())
//    YMIR_TEST(test_model_gene_usage())
//    YMIR_TEST(test_model_vj_maag())
//    YMIR_TEST(test_model_vdj_maag())
    

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