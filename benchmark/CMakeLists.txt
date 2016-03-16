add_executable (Benchmark benchmark.cpp)
add_executable (Bench_Inference bench_inference.cpp)

include_directories(../Ymir)
include_directories(../include)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build")