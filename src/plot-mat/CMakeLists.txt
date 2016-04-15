project(e-scatter/plot-mat)
add_executable(plot-mat main.cc
	../cdsem/material.cc)
target_link_libraries(plot-mat common)
