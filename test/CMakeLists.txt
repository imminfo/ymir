option(TestFolderPath "/Users/vdn/Projects/ymir/test/data/" ON)


include_directories(../Ymir)
include_directories(../include)
include_directories(/usr/local/include/)

add_executable(Test_ParamVec_and_VDJGenes test_paramvec_genes.cpp)
add_executable(Test_Alignment_Clone test_alignment_clone.cpp)
add_executable(Test_Aligner test_aligner.cpp)
add_executable(Test_IO test_io.cpp)
add_executable(Test_MarkovChain test_markov_chain.cpp)
add_executable(Test_MAAG test_maag.cpp)
add_executable(Test_Model test_model.cpp)


enable_testing()

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build")

add_test(MyTest Ymir-test)