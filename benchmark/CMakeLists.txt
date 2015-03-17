add_executable (Ymir-benchmark benchmark.cpp)

include_directories(../src)
include_directories(../include)

enable_testing()

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build")