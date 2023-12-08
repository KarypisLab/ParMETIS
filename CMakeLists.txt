cmake_minimum_required(VERSION 3.14...3.28)
project(ParMETIS LANGUAGES C)

# --- need METIS first
include(FetchContent)

FetchContent_Declare(GKlib
GIT_REPOSITORY https://github.com/scivision/GKlib.git
GIT_TAG 44630ca
GIT_SHALLOW true
)
FetchContent_MakeAvailable(GKlib)

FetchContent_Declare(metis_fc
GIT_REPOSITORY https://github.com/scivision/metis.git
GIT_TAG fd7799b
GIT_SHALLOW true
)

FetchContent_MakeAvailable(metis_fc)

if(NOT DEFINED REALTYPEWIDTH)
  set(REALTYPEWIDTH 32)
endif()
if(NOT DEFINED IDXTYPEWIDTH)
  set(IDXTYPEWIDTH 32)
endif()
add_compile_definitions(REALTYPEWIDTH=${REALTYPEWIDTH} IDXTYPEWIDTH=${IDXTYPEWIDTH})

message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION} CMake ${CMAKE_VERSION}
REALTYPEWIDTH=${REALTYPEWIDTH} IDXTYPEWIDTH=${IDXTYPEWIDTH}")
message(STATUS "METIS ${metis_fc_SOURCE_DIR}")
message(STATUS "GKlib ${gklib_SOURCE_DIR}")

find_package(MPI REQUIRED)

include(./conf/gkbuild.cmake)

install(FILES include/parmetis.h TYPE INCLUDE)

add_subdirectory(libparmetis)
add_subdirectory(programs)

# This is for testing during development and is not being distributed
#add_subdirectory(test)
