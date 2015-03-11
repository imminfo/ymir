add_executable (Ymir-test test.cpp)

include_directories(../src)
include_directories(../include)

enable_testing()

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build")

add_test(MyTest Ymir-test)