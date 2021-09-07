cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0048 NEW)

project(mmtdt
        VERSION 0.0.0
        DESCRIPTION "Mass mapping with transdimensional trees"
        LANGUAGES CXX C)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(PROJECT_LIB lib${PROJECT_NAME})

set(CMAKE_VERBOSE_MAKEFILE OFF)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option(tests "Enable testing" on)
option(install "Install targets on system" on)

include(ExternalProject)
include(FetchContent)
include(GetTDTBase)

set(SOURCES src/mmobservations.cpp
            src/proposals.cpp
            src/globalprop.cpp
            src/birthprop.cpp
            src/deathprop.cpp
            src/valueprop.cpp
            src/rng.cpp
            src/hierarchicalmodel.cpp
            src/logging.cpp
            src/utils.cpp)

add_library(${PROJECT_LIB} STATIC ${SOURCES})
set_target_properties(${PROJECT_LIB} PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
add_dependencies(${PROJECT_LIB} ${TDTBASE})

find_package(GSL REQUIRED)
find_package(GMP REQUIRED)
find_package(FFTW3 REQUIRED)
target_link_libraries(${PROJECT_LIB}
                      GSL::gsl
                      ${GMP_LIBRARY}
                      ${TDTBASE_LIBS}
                      FFTW3::FFTW3)

include_directories(include ${TDTBASE_INCLUDE}
                            ${GMP_INCLUDE_DIR})

if(tests)
    enable_testing()
    include(CTest)
    FetchContent_Declare(googletest
                        GIT_REPOSITORY https://github.com/google/googletest.git
                        GIT_TAG master
                        INSTALL_COMMAND "")
    FetchContent_MakeAvailable(googletest)
    add_subdirectory(tests)
    configure_file(data/Bolshoi_7_clean_256.txt ${CMAKE_CURRENT_BINARY_DIR}/Bolshoi_7_clean_256.txt)
    configure_file(data/tutorial_prior.txt ${CMAKE_CURRENT_BINARY_DIR}/tutorial_prior.txt)
endif()

add_executable(${PROJECT_NAME} src/massmapping.cpp)
target_link_libraries(${PROJECT_NAME} ${PROJECT_LIB})

add_executable(${PROJECT_NAME}_pp src/postprocess.cpp)
target_link_libraries(${PROJECT_NAME}_pp ${PROJECT_LIB})

if(install)
    install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_pp RUNTIME)
    install(TARGETS ${PROJECT_LIB} LIBRARY)
endif()