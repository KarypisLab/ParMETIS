include_directories(.)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-error")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-error")

# Create programs.
add_executable(pm_ptest ptest.c io.c adaptgraph.c)
add_executable(pm_mtest mtest.c io.c)
add_executable(pm_parmetis parmetis.c io.c adaptgraph.c)
add_executable(pm_pometis pometis.c io.c)
add_executable(pm_dglpart dglpart.c)
add_executable(pm_dglpart2 dglpart2.c)
add_executable(pm_dglpart3 dglpart3.c)

# Link with the required libraries 
foreach(prog pm_ptest pm_mtest pm_parmetis pm_pometis pm_dglpart pm_dglpart2 pm_dglpart3)
  target_link_libraries(${prog} parmetis metis GKlib m)
endforeach(prog)

install(TARGETS pm_ptest pm_mtest pm_parmetis pm_pometis pm_dglpart pm_dglpart2 pm_dglpart3
  RUNTIME DESTINATION bin)
