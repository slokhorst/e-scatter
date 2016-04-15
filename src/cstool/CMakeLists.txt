project(e-scatter/cstool)
add_executable(cstool main.cc cross-section.cc 
	../cdsem/material.cc
	../compile-mat/parser.cc)
target_link_libraries(cstool common muparser)
