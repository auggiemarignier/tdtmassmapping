cmake_minimum_required(VERSION 3.18)

set(LIBS gtest_main ${PROJECT_LIB})

add_executable(test_proposals test_proposals.cpp)
add_executable(test_mmobservations test_mmobservations.cpp)

target_link_libraries(test_proposals PUBLIC ${LIBS})
target_link_libraries(test_mmobservations PUBLIC ${LIBS})

add_test(NAME test_proposals COMMAND test_proposals)
add_test(NAME test_mmobservations COMMAND test_mmobservations)

configure_file(../data/tutorial_prior.txt ${CMAKE_CURRENT_BINARY_DIR}/tutorial_prior.txt)
