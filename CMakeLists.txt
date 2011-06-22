cmake_minimum_required(VERSION 2.6)
#project(optpp)
project(optpp C CXX Fortran)

include(FortranCInterface)

find_library(npsol_lib libnpsol.a /usr/local/lib/libnpsol.a)
include_directories(
  ${optpp_SOURCE_DIR}/include
  ${optpp_SOURCE_DIR}/newmat11)

find_package(BLAS)
if(BLAS_FOUND)
  set(HAVE_BLAS TRUE)
endif()

if(NOT BUILT_BY_ACRO)
  option(USE_MPI "Should OPTPP use MPI?" OFF)
  if(USE_MPI)
    find_package(MPI)
  endif()
else()
  set(DAKOTA_OPTPP TRUE)
endif()
if(Dakota_SOURCE_DIR)
  set(DAKOTA_OPTPP TRUE)
endif()

if(MPI_FOUND)
  set(WITH_MPI TRUE)
  include_directories(${MPI_INCLUDE_PATH})
endif()

option(OPTPP_SHARED "Turn ON if you have a shared file system" OFF)
if(OPTPP_SHARED)
  set(SHARED TRUE)
  if(NOT MPI_FOUND)
    message(FATAL_ERROR " optpp cannot be built on a shared file system without MPI.  Please set WITH_MPI to ON")
  endif()
endif()

include(CheckFunctionExists)
CHECK_FUNCTION_EXISTS(gettimeofday HAVE_GETTIMEOFDAY)
CHECK_FUNCTION_EXISTS(times HAVE_TIMES)

include(CheckIncludeFiles)
CHECK_INCLUDE_FILES(sys/param.h HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILES(sys/resource.h HAVE_SYS_RESOURCE_H)
CHECK_INCLUDE_FILES(sys/times.h HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILES(sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
# WJB - ToDo: ask Zack to review this change a recommend a better approach
#             to avoid macro redefinition warnings (see OPT++_config.h too).
if(HAVE_SYS_TYPES_H)
  add_definitions("-DHAVE_SYS_TYPES_H")
endif(HAVE_SYS_TYPES_H)


try_compile(HAVE_NAMESPACES
  ${optpp_BINARY_DIR}/config/have_namespaces
  ${optpp_SOURCE_DIR}/config/have_namespaces.cxx
)
try_compile(HAVE_STD
  ${optpp_BINARY_DIR}/config/have_std
  ${optpp_SOURCE_DIR}/config/have_std.cxx
)

configure_file(
  ${optpp_SOURCE_DIR}/config/OPT++_config_cmake.h.in
  ${optpp_BINARY_DIR}/config/OPT++_config.h
)
include_directories(${optpp_BINARY_DIR}/config)
add_definitions(-DHAVE_CONFIG_H)

add_subdirectory(newmat11)
add_subdirectory(src)
enable_testing()
add_subdirectory(tests)
