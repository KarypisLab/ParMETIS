# Find sources.
set(parmetis_sources diffutil.c	initmsection.c	move.c		rename.h	wave.c
akwayfm.c	frename.c	initpart.c	msetup.c	renumber.c	weird.c
ametis.c	gklib.c		kmetis.c	node_refine.c	rmetis.c	wspace.c
balancemylink.c	gklib_defs.h	kwayrefine.c	ometis.c	selectq.c	xyzpart.c
comm.c		gklib_rename.h	macros.h	parmetislib.h	serial.c
csrmatch.c	gkmetis.c	match.c		proto.h		stat.c
ctrl.c		gkmpi.c		mdiffusion.c	pspases.c	struct.h
debug.c		graph.c		mesh.c		redomylink.c	timer.c
defs.h		initbalance.c	mmetis.c	remap.c		util.c)

# Create libparmetis
add_library(parmetis ${parmetis_sources})
target_include_directories(parmetis PUBLIC
"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR};${PROJECT_SOURCE_DIR}/include>"
$<INSTALL_INTERFACE:include>
PRIVATE
${gklib_SOURCE_DIR}
${metis_fc_SOURCE_DIR}/include
)
target_link_libraries(parmetis PRIVATE MPI::MPI_C)

if(SHARED)
  target_link_libraries(parmetis PRIVATE metis GKlib)
endif()

install(TARGETS parmetis)
