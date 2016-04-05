project(e-scatter/compile-mat)
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})
add_executable(compile-mat main.cc 
	../common/spline.cc ../common/archive.cc ../common/polynomial.cc ../common/interpolate.cc ../common/optional.cc 
	../cdsem/material.cc 
	../cstool/cross-section.cc 
	../cpl/xml.cc ../cpl/parser.cc ../cpl/text.cc ../cpl/mathexpr.cc)
target_link_libraries(compile-mat ${Boost_LIBRARIES})
