cmake_minimum_required(VERSION 2.8.4)

project(Ymir)

#add_definitions(-DUSE_OMP)

#if (USE_OMP)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -fopenmp")

# /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++
# /usr/local/Cellar/clang-omp/2015-04-01/bin/clang-omp++

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=gnu++11")
endif()

include_directories(include)
include_directories(/usr/local/include/)

add_subdirectory(test)
add_subdirectory(benchmark)
add_subdirectory(scripts)