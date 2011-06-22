add_library(rosen rosen.C)
target_link_libraries(rosen ${test_libs})

add_library(hockfcns hockfcns.C)
target_link_libraries(hockfcns ${test_libs})
