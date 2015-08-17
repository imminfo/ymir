add_executable (Test test.cpp)

include_directories(../Ymir)
include_directories(../include)
include_directories(/usr/local/include/)

enable_testing()

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build")

add_test(MyTest Ymir-test)