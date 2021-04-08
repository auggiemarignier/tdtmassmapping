cmake_minimum_required(VERSION 3.18)

set(SOURCES tests.cpp)

add_executable(tests ${SOURCES})
target_link_libraries(tests PUBLIC gtest_main ${PROJECT_LIB})

add_test(NAME tests COMMAND tests)