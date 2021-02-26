cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0048 NEW)
set(CMAKE_VERBOSE_MAKEFILE ON)

project(mmtdt
        VERSION 0.0.0
        DESCRIPTION "Mass mapping with transdimensional trees"
        LANGUAGES CXX)

include(ExternalProject)

set(TDTBASE tdtbase)
ExternalProject_Add(${TDTBASE}
    GIT_REPOSITORY git@github.com:rhyshawkins/TDTbase.git
    GIT_SHALLOW ON
    GIT_SUBMODULES ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE ON
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(${TDTBASE} SOURCE_DIR)