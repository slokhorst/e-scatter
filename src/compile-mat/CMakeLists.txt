project(e-scatter/compile-mat)
add_executable(compile-mat main.cc 
	../common/spline.cc ../common/archive.cc ../common/polynomial.cc ../common/interpolate.cc ../common/optional.cc 
	../cdsem/material.cc 
	../cstool/cross-section.cc 
	../cpl/xml.cc ../cpl/parser.cc ../cpl/text.cc ../cpl/mathexpr.cc)
target_link_libraries(compile-mat boost_program_options)
