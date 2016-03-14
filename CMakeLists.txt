cmake_minimum_required(VERSION 2.8.4)

project(Ymir)

#set(CMAKE_C_COMPILER /Users/vdn/Downloads/clang-llvm/bin/clang)
#set(CMAKE_CXX_COMPILER /Users/vdn/Downloads/clang-llvm/bin/clang++)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
endif()

#set(SOURCE_FILES Ymir/src/main.cpp)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build")
#include(CTest)

include_directories(include)
include_directories(/usr/local/include/)

add_subdirectory(test)
add_subdirectory(benchmark)
#add_subdirectory(scripts)