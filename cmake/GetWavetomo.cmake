set(TDTWAVETOMO wavetomo)
ExternalProject_Add(${TDTWAVETOMO}
    GIT_REPOSITORY git@github.com:rhyshawkins/TDTWavetomo2D.git
    GIT_SHALLOW ON
    GIT_SUBMODULES ""
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make TDTBASE=${TDTBASE_SRC} TRAVELTIMEBASE=${TRAVELTIME_SRC} CXX=${CMAKE_CXX_COMPILER}
    BUILD_IN_SOURCE ON
    INSTALL_COMMAND ""
    DEPENDS ${TRAVELTIME}
)
ExternalProject_Get_Property(${TDTWAVETOMO} SOURCE_DIR)
set(TDTWAVETOMO_INCLUDE ${SOURCE_DIR})

set(TDTWAVETOMO_LIB ${TDTWAVETOMO}_LIB)
add_library(${TDTWAVETOMO_LIB} OBJECT IMPORTED GLOBAL)
set_property(TARGET ${TDTWAVETOMO_LIB} PROPERTY IMPORTED_OBJECTS
            ${SOURCE_DIR}/birth.o
            ${SOURCE_DIR}/birthslice.o
            ${SOURCE_DIR}/death.o
            ${SOURCE_DIR}/deathslice.o
            ${SOURCE_DIR}/global.o
            ${SOURCE_DIR}/globalslice.o
            ${SOURCE_DIR}/hierarchical.o
            ${SOURCE_DIR}/hierarchicalmodel.o
            ${SOURCE_DIR}/hierarchicalprior.o
            ${SOURCE_DIR}/hierarchicalpriorslice.o
            ${SOURCE_DIR}/hierarchicalslice.o
            ${SOURCE_DIR}/ptexchange.o
            ${SOURCE_DIR}/ptexchangeslice.o
            ${SOURCE_DIR}/resample.o
            ${SOURCE_DIR}/rng.o
            ${SOURCE_DIR}/value.o
            ${SOURCE_DIR}/valueslice.o
            ${SOURCE_DIR}/volume.o
            ${SOURCE_DIR}/wavetomo2dexception.o
            ${SOURCE_DIR}/wavetomo2dobservations.o
            ${SOURCE_DIR}/wavetomo2dutil.o
            )