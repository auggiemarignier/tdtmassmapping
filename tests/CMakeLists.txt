cmake_minimum_required(VERSION 3.18)

set(LIBS gtest_main ${PROJECT_LIB})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

foreach(test test_proposals
             test_mmobservations
             test_globalprop
             test_utils
             test_ggd_prior)
    add_executable(${test} ${test}.cpp)
    target_link_libraries(${test} PUBLIC ${LIBS})
    add_test(NAME ${test} COMMAND ${test})
endforeach()

