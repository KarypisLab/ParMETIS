# Include directories for library code.
include_directories(.)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")

# Find sources.
file(GLOB parmetis_sources *.c)

# Create libparmetis
add_library(parmetis ${ParMETIS_LIBRARY_TYPE} ${parmetis_sources})

install(TARGETS parmetis
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib)
