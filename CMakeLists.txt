cmake_minimum_required(VERSION 2.8.4)
project(Ymir)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Xlang -fopenmp=libiomp5")

set(SOURCE_FILES src/main.cpp)
add_executable(Ymir ${SOURCE_FILES})

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build")

include_directories(include)
add_subdirectory(test)
add_subdirectory(benchmark)