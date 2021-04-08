cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0048 NEW)

project(mmtdt
        VERSION 0.0.0
        DESCRIPTION "Mass mapping with transdimensional trees"
        LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option(tests "Enable testing" on)

include(ExternalProject)
include(FetchContent)
include(dependencies)

set(SOURCES src/massmapping.cpp
            src/proposals.cpp)

add_executable(mmtdt ${SOURCES}
                     $<TARGET_OBJECTS:${TDTWAVETOMO_LIB}>
                     ${TRAVELTIMEEXCEPTION})
set_property(TARGET mmtdt PROPERTY CXX_STANDARD 11)
add_dependencies(mmtdt ${TDTWAVETOMO})

find_package(GSL REQUIRED)
find_package(GMP REQUIRED)
target_link_libraries(mmtdt
                    GSL::gsl
                    ${GMP_LIBRARY}
                    ${TRAVELTIME_SRC}/libtraveltime2d.a
                    ${TDTBASE_LIBS})

include_directories(include ${TDTWAVETOMO_INCLUDE}
                            ${TDTBASE_INCLUDE}
                            ${TRAVELTIME_INCLUDE})

if(tests)
    enable_testing()
    FetchContent_Declare(googletest
                        GIT_REPOSITORY git@github.com:google/googletest.git
                        GIT_TAG master)
    FetchContent_MakeAvailable(googletest)
endif()
