project(e-scatter/plot-mat)
add_executable(plot-mat main.cc
	../common/archive.cc ../common/spline.cc ../common/polynomial.cc ../common/interpolate.cc ../common/optional.cc
	../cdsem/material.cc
	../cpl/text.cc ../cpl/parser.cc)
target_link_libraries(plot-mat)
