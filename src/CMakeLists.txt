if(WIN32)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)

add_subdirectory(Base)
add_subdirectory(Constraints)
add_subdirectory(GSS)
add_subdirectory(Newton)
add_subdirectory(PDS)
add_subdirectory(Utils)

add_library(optpp globals.C)
target_link_libraries(optpp PUBLIC base constraints gss newton pds utils)
target_compile_definitions(optpp PUBLIC HAVE_OPTPP_CONFIG_H)
set_target_properties(optpp PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${optpp_BINARY_DIR}/lib
  ARCHIVE_OUTPUT_DIRECTORY ${optpp_BINARY_DIR}/lib)
  
install(TARGETS optpp
    EXPORT optpp-config
	DESTINATION lib
	INCLUDES DESTINATION include/optpp
)
