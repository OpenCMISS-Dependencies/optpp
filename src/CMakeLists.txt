add_subdirectory(Base)
add_subdirectory(Constraints)
add_subdirectory(GSS)
add_subdirectory(Newton)
add_subdirectory(PDS)
add_subdirectory(Utils)

add_library(opt globals.C)
target_link_libraries(opt base constraints gss newton pds utils)
set_target_properties(opt PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${optpp_BINARY_DIR}/lib
  ARCHIVE_OUTPUT_DIRECTORY ${optpp_BINARY_DIR}/lib)
