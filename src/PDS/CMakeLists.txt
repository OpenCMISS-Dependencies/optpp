set(libpds_srcs
  create_scheme.C dqrdc.c make_search.c pds.C pdschk.C pdscld.c pdscom.c
  pdsdgn.c pdsdone.c pdseql.c pdsget.c pdshrk.c pdsinit.C pdslen.c pdsopt.C
  pdsrgt.c pdsstep.C pdsupd.c pdswork.C sort.c writes.c
)

if(MPI_FOUND)
  list(APPEND libpds_srcs pdsglb.c pdsgop.c pdswap.c) 
endif()

add_library(pds ${libpds_srcs})

set(pds_libs utils)
if(BLAS_FOUND)
  list(APPEND pds_libs ${BLAS_LIBRARIES})
endif(BLAS_FOUND)

target_link_libraries(pds ${pds_libs})

INSTALL(TARGETS pds
	ARCHIVE DESTINATION lib/optpp
)
