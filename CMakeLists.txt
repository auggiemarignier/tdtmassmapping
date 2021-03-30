cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0048 NEW)

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
ExternalProject_Get_property(${TDTBASE} SOURCE_DIR)
set(TDTBASESRC ${SOURCE_DIR})

set(TRAVELTIME traveltime2d)
ExternalProject_Add(${TRAVELTIME}
    GIT_REPOSITORY git@github.com:rhyshawkins/traveltime2d.git
    GIT_SHALLOW ON
    GIT_SUBMODULES ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_IN_SOURCE ON
    INSTALL_COMMAND ""
    DEPENDS ${TDTBASE}
)
ExternalProject_Get_property(${TRAVELTIME} SOURCE_DIR)
set(TRAVELTIMESRC ${SOURCE_DIR}/..)

set(TDTWAVETOMO wavetomo)
ExternalProject_Add(${TDTWAVETOMO}
    GIT_REPOSITORY git@github.com:rhyshawkins/TDTWavetomo2D.git
    GIT_SHALLOW ON
    GIT_SUBMODULES ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make TDTBASE=${TDTBASESRC} TRAVELTIMEBASE=${TRAVELTIMESRC} CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE ON
    INSTALL_COMMAND ""
    DEPENDS ${TRAVELTIME}
)

include_directories(include)
set(SOURCES src/massmapping.cpp)

add_executable(mmtdt ${SOURCES})