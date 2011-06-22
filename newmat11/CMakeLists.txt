add_library(newmat
  fft.C evalue.C submat.C cholesky.C hholder.C sort.C newmatrm.C jacobi.C svd.C
  newmat8.C newmat7.C newmat6.C newmat5.C newmat3.C newmat4.C newmat2.C
  newmat1.C bandmat.C myexcept.C newmatex.C newfft.C
)

set_target_properties(newmat PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${optpp_BINARY_DIR}/lib
  ARCHIVE_OUTPUT_DIRECTORY ${optpp_BINARY_DIR}/lib)
  
INSTALL(TARGETS newmat
	ARCHIVE DESTINATION lib/optpp
)

# install the header files
SET(NEWMAT_INCLUDES controlw.h    include.h       myexcept.h
                  newmatap.h    newmat.h        newmatio.h
                  newmatnl.h    newmatrc.h      newmatrm.h
                  precisio.h    solution.h      tmt.h)

INSTALL(FILES ${NEWMAT_INCLUDES}
	DESTINATION include/optpp
)
