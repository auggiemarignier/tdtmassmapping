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
set(TDTBASE_SRC ${SOURCE_DIR})
list(APPEND TDTBASE_INCLUDE
    ${TDTBASE_SRC}/hnk
    ${TDTBASE_SRC}/log
    ${TDTBASE_SRC}/oset
    ${TDTBASE_SRC}/sphericalwavelet
    ${TDTBASE_SRC}/tracking
    ${TDTBASE_SRC}/wavelet
    ${TDTBASE_SRC}/wavetree
    )