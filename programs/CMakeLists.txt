include_directories(.)

# Create programs.
add_executable(pm_ptest ptest.c io.c adaptgraph.c)
add_executable(pm_mtest mtest.c io.c)
add_executable(pm_parmetis parmetis.c io.c adaptgraph.c)
add_executable(pm_pometis pometis.c io.c)
add_executable(pm_dglpart dglpart.c)

# Link with the required libraries 
foreach(prog pm_ptest pm_mtest pm_parmetis pm_pometis pm_dglpart)
  target_link_libraries(${prog} parmetis metis GKlib m)
endforeach(prog)

install(TARGETS pm_ptest pm_mtest pm_parmetis pm_pometis pm_dglpart
  RUNTIME DESTINATION bin)
