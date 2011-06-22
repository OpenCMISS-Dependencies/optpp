add_library(base 
  Appl_Data.C Appl_Data_NPSOL.C backtrack.C dogleg.C FDNLF1.C linesearch.C
  LSQNLF.C mcsrch.C NLF0.C NLF1.C NLF2.C NLP0.C NLP1.C NLP2.C NLP.C TOLS.C
  trustpds.C trustregion.C
)

set(base_libs constraints pds)
if(MPI_FOUND)
  list(APPEND base_libs ${MPI_LIBRARY})
endif()
target_link_libraries(base ${base_libs})

INSTALL(TARGETS base
	ARCHIVE DESTINATION lib/optpp
)
