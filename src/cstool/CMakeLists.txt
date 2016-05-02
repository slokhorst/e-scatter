project(e-scatter/cstool)
include(FindPkgConfig)

pkg_search_module(muParser muparser)
if(NOT muParser_FOUND)
	message(WARNING "missing libraries, not building ${PROJECT_NAME}")
	return()
endif(NOT muParser_FOUND)

include_directories(${muParser_INCLUDE_DIRS})
add_executable(cstool main.cc cross-section.cc 
	../cdsem/material.cc
	../compile-mat/parser.cc)
target_link_libraries(cstool common ${muParser_LIBRARIES})
