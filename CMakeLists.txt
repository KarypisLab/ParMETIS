cmake_minimum_required(VERSION 3.21...4.3)

if(NOT CMAKE_BUILD_TYPE AND NOT DEFINED ENV{CMAKE_BUILD_TYPE})
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release")
endif()

project(PARMETIS LANGUAGES C VERSION 1.0.0)

# --- need METIS first
include(GNUInstallDirs)
include(FetchContent)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND PARMETIS_IS_TOP_LEVEL)
  set(CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/local" CACHE PATH "install prefix" FORCE)
endif()

message(STATUS "${PROJECT_NAME} ${PROJECT_VERSION} CMake ${CMAKE_VERSION} install prefix ${CMAKE_INSTALL_PREFIX}")

set(metis_url "https://github.com/scivision/METIS/archive/456244e99d7f0027706033d18bc1b3424d3507cc.tar.gz")
FetchContent_Declare(METIS URL ${metis_url})
FetchContent_MakeAvailable(METIS)

set(CMAKE_C_STANDARD 99)

add_compile_definitions(
"REALTYPEWIDTH=$<IF:$<BOOL:${REALTYPEWIDTH}>,${REALTYPEWIDTH},32>"
"IDXTYPEWIDTH=$<IF:$<BOOL:${IDXTYPEWIDTH}>,${IDXTYPEWIDTH},32>"
)

find_package(MPI REQUIRED)

include(./conf/gkbuild.cmake)

install(FILES include/parmetis.h TYPE INCLUDE)

add_library(parmetis)
target_include_directories(parmetis PUBLIC
"$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include;${PROJECT_SOURCE_DIR}/libparmetis;${metis_BINARY_DIR}>"
$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(parmetis PRIVATE MPI::MPI_C METIS::METIS)

if(SHARED)
  target_link_libraries(parmetis PRIVATE metis GKlib)
endif()

add_library(PARMETIS::PARMETIS ALIAS parmetis)

install(TARGETS parmetis EXPORT ${PROJECT_NAME}-targets)

add_subdirectory(libparmetis)
add_subdirectory(programs)

# This is for testing during development and is not being distributed
#add_subdirectory(test)

include(cmake/install.cmake)

file(GENERATE OUTPUT .gitignore CONTENT "*")
