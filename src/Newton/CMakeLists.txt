set(libnewton_srcs
  OptBaNewton.C OptBaQNewton.C OptBCEllipsoid.C OptBCFDNewton.C OptBCNewton.C
  OptBCNewtonLike.C OptBCQNewton.C OptCG.C OptConstrFDNewton.C OptConstrNewton.C
  OptConstrNewtonLike.C OptConstrQNewton.C OptDHNIPS.C OptFDNewton.C OptFDNIPS.C
  OptLBFGS.C OptNewton.C OptNewtonLike.C OptNIPS.C OptNIPSLike.C OptQNewton.C
  OptQNIPS.C ../globals.C
)

if(npsol_lib)
  list(APPEND libnewton_srcs OptNPSOL.C npsol_setup.c)
endif()

add_library(newton ${libnewton_srcs})
target_link_libraries(newton newmat base)

INSTALL(TARGETS newton
    EXPORT optpp-config
	ARCHIVE DESTINATION lib
	INCLUDES DESTINATION include/optpp
)