cmake_minimum_required(VERSION 3.18)

set(LIBS gtest_main ${PROJECT_LIB})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)

add_executable(test_proposals test_proposals.cpp)
add_executable(test_mmobservations test_mmobservations.cpp)
add_executable(test_globalprop test_globalprop.cpp)

target_link_libraries(test_proposals PUBLIC ${LIBS})
target_link_libraries(test_mmobservations PUBLIC ${LIBS})
target_link_libraries(test_globalprop PUBLIC ${LIBS})

add_test(NAME test_proposals COMMAND test_proposals)
add_test(NAME test_mmobservations COMMAND test_mmobservations)
add_test(NAME test_globalprop COMMAND test_globalprop)
