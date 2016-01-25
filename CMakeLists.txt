cmake_minimum_required(VERSION 3.0)
project(optpp VERSION 681 LANGUAGES C CXX)
enable_language(Fortran OPTIONAL)
set(CMAKE_POSITION_INDEPENDENT_CODE YES)

option(USE_NPSOL "Do you want to use the npsol library?" OFF)
if(USE_NPSOL)
  find_library(npsol_lib libnpsol.a /usr/local/lib/libnpsol.a)
endif(USE_NPSOL)

include_directories(
  ${optpp_SOURCE_DIR}/include
  ${optpp_SOURCE_DIR}/newmat11)

option(USE_EXTERNAL_BLAS "Allow to use external BLAS libs?" OFF)
if(USE_EXTERNAL_BLAS)
  include(FortranCInterface)
  find_package(BLAS ${BLAS_VERSION} QUIET)
  if(BLAS_FOUND)
    set(HAVE_BLAS TRUE)
  endif()
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
add_definitions(-DHAVE_OPTPP_CONFIG_H)

INSTALL(FILES ${CMAKE_BINARY_DIR}/config/OPT++_config.h
	DESTINATION include/optpp
)

SET(CMAKE_DEBUG_POSTFIX "d")

add_subdirectory(newmat11)
add_subdirectory(src)
#enable_testing()
#add_subdirectory(tests)

# install the header files
SET(OPT_INCLUDES
  include/abbrev_copyright.h    include/AppLauncher.h       
  include/Appl_Data.h           include/Appl_Data_NPSOL.h   
  include/BoolVector.h          include/BoundConstraint.h   
  include/cblas.h               include/CGProblem.h         
  include/common.h              include/CompoundConstraint.h
  include/ConstraintBase.h      include/Constraint.h        
  include/GenSetBase.h          include/GenSetBox2d.h       
  include/GenSet.h              include/GenSetMin.h         
  include/GenSetStd.h           include/globals.h           
  include/ioformat.h            include/LinearConstraint.h  
  include/LinearEquation.h      include/LinearInequality.h  
  include/LSQNLF.h              include/NewtonProblem.h     
  include/NIPSProblem.h         include/NLFAPP.h            
  include/NLF.h                 include/NLP0.h              
  include/NLP1.h                include/NLP2.h              
  include/NLPBase.h             include/NLP.h               
  include/NonLinearConstraint.h include/NonLinearEquation.h 
  include/NonLinearInequality.h include/NPSOLProblem.h      
  include/OptBaNewton.h         include/OptBaQNewton.h      
  include/OptBCEllipsoid.h      include/OptBCFDNewton.h     
  include/OptBCNewton.h         include/OptBCNewtonLike.h   
  include/OptBCQNewton.h        include/OptCG.h             
  include/OptConstrFDNewton.h   include/OptConstrNewton.h   
  include/OptConstrNewtonLike.h include/OptConstrQNewton.h  
  include/OptDHNIPS.h           include/OptDirect.h         
  include/OptFDNewton.h         include/OptFDNIPS.h         
  include/OptGSS.h              include/Opt.h               
  include/OptLBFGS.h            include/OptNewton.h         
  include/OptNewtonLike.h       include/OptNIPS.h           
  include/OptNIPSLike.h         include/OptNPSOL.h          
  include/Opt_PARAMS.h          include/OptPDS.h            
  include/OptppArray.h          include/OptppExceptions.h   
  include/OptppFatalError.h     include/OptppSmartPtr.h     
  include/OptQNewton.h          include/OptQNIPS.h          
  include/pds.h                 include/PDSProblem.h        
  include/Problem.h             include/proto.h             
  include/TOLS.h                include/VariableList.h
)

INSTALL(FILES ${OPT_INCLUDES}
	DESTINATION include/optpp
)

install(EXPORT optpp-config
	DESTINATION ${PACKAGE_CONFIG_DIR})
include(CMakePackageConfigHelpers)
WRITE_BASIC_PACKAGE_VERSION_FILE(${CMAKE_CURRENT_BINARY_DIR}/optpp-config-version.cmake
    COMPATIBILITY AnyNewerVersion)
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/optpp-config-dependencies.cmake 
    "include(CMakeFindDependencyMacro)\r\n"
	"set(${PROJECT_NAME}_IMPORT_PREFIX \${_IMPORT_PREFIX})\r\n"
)
if (USE_EXTERNAL_BLAS AND HAVE_BLAS)
    file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/optpp-config-dependencies.cmake "find_dependency(BLAS)\r\n")
endif()
file(APPEND ${CMAKE_CURRENT_BINARY_DIR}/optpp-config-dependencies.cmake "set(_IMPORT_PREFIX \${${PROJECT_NAME}_IMPORT_PREFIX})")    
install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/optpp-config-version.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/optpp-config-dependencies.cmake
    DESTINATION ${PACKAGE_CONFIG_DIR})