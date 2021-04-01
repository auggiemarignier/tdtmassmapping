cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0048 NEW)

project(mmtdt
        VERSION 0.0.0
        DESCRIPTION "Mass mapping with transdimensional trees"
        LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

include(ExternalProject)
include(dependencies)

set(SOURCES src/massmapping.cpp)

add_executable(mmtdt ${SOURCES})
set_property(TARGET mmtdt PROPERTY CXX_STANDARD 11)
add_dependencies(mmtdt ${TDTWAVETOMO})

find_package(GSL REQUIRED)
target_link_libraries(mmtdt GSL::gsl)

include_directories(include ${TDTWAVETOMO_INCLUDE}
                            ${TDTBASE_INCLUDE}
                            ${TRAVELTIME_INCLUDE})