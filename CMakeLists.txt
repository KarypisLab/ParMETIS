cmake_minimum_required(VERSION 3.10)
project(ParMETIS C)


# Search for MPI.
# GK commented this out as it seems to be creating problems
# include(FindMPI)
# if(NOT MPI_FOUND)
#   message(FATAL_ERROR "mpi is not found")
# endif()
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")


# Prepare libraries.
if(SHARED)
  set(ParMETIS_LIBRARY_TYPE SHARED)
else()
  set(ParMETIS_LIBRARY_TYPE STATIC)
endif()

include(./conf/gkbuild.cmake)

# List of paths that the compiler will search for header files.
# i.e., the -I equivalent
include_directories(include)
include_directories(${MPI_INCLUDE_PATH})
include_directories(${GKLIB_PATH}/include)
include_directories(${METIS_PATH}/include)
include_directories(${CMAKE_INSTALL_PREFIX}/include)

# List of paths that the compiler will search for library files.
# i.e., the -L equivalent
link_directories(${GKLIB_PATH}/lib)
link_directories(${METIS_PATH}/lib)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

# List of directories that cmake will look for CMakeLists.txt
add_subdirectory(include)
add_subdirectory(libparmetis)
add_subdirectory(programs)

# This is for testing during development and is not being distributed
#add_subdirectory(test)
