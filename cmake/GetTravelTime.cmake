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
set(TRAVELTIME_SRC ${SOURCE_DIR}/..)
set(TRAVELTIME_INCLUDE ${SOURCE_DIR})