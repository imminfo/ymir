cmake_minimum_required(VERSION 2.8.4)

project(Ymir)


set(USE_OMP NO CACHE BOOL "Use OpenMP or not.")

if (USE_OMP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()


set(OMP_THREADS -1 CACHE STRING "Number of OpenMP threads to use.")


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
add_subdirectory(scripts)