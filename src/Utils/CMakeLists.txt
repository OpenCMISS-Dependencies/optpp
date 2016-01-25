set(libutils_srcs
  BoolVector.C file_cutils.c ioformat.C mcholesky.C OptppExceptions.C
  OptppFatalError.C print.C timers.c
)

if(NOT USE_EXTERNAL_BLAS OR NOT HAVE_BLAS)
  list(APPEND libutils_srcs linalg.c)
endif()

INCLUDE(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(times HAVE_TIMES)
if(HAVE_TIMES)
  add_definitions(-DHAVE_TIMES)
endif()
INCLUDE(CheckIncludeFiles)
CHECK_INCLUDE_FILES("sys/times.h" HAVE_SYS_TIMES_H)
if(HAVE_SYS_TIMES_H)
  add_definitions(-DHAVE_SYS_TIMES_H)
endif()
CHECK_INCLUDE_FILES("sys/time.h" HAVE_SYS_TIME_H)
if(HAVE_SYS_TIME_H)
  add_definitions(-DHAVE_SYS_TIME_H)
endif()

add_library(utils ${libutils_srcs})

INSTALL(TARGETS utils
    EXPORT optpp-config
	DESTINATION lib
	INCLUDES DESTINATION include/optpp
)
