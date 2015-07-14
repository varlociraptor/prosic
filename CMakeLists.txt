cmake_minimum_required (VERSION 2.6)
project (sm_caller)

# The version number: 
set (sm_caller_VERSION_MAJOR 1)
set (sm_caller_VERSION_MINOR 0)

file(GLOB SOURCES src/*.c)

# Set the output directory for the build executables
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_BINARY_DIR}/config.h"
)

# add the binary tree to the search path for include files
# so that we will find config.h
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/src")

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}) 

FIND_PACKAGE(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GSL_LIBRARIES})

FIND_PACKAGE(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIR})
set(LIBS ${LIBS} ${GMP_LIBRARIES})


add_executable(sm_caller ${SOURCES})
target_link_libraries (sm_caller ${LIBS})

install (TARGETS sm_caller DESTINATION "${PROJECT_SOURCE_DIR}/bin")
install (FILES "${PROJECT_BINARY_DIR}/config.h" DESTINATION "${PROJECT_SOURCE_DIR}/include")
install (FILES "src/smc_input.h" DESTINATION "${PROJECT_SOURCE_DIR}/include")
install (FILES "src/smc_likelihood.h" DESTINATION "${PROJECT_SOURCE_DIR}/include")
