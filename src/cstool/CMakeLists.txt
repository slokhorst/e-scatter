project(e-scatter/cstool)
find_package(PkgConfig REQUIRED)
set (CMAKE_CXX_STANDARD 14)

pkg_search_module(muParser muparser)
if(NOT muParser_FOUND)
	message(WARNING "missing libraries, not building ${PROJECT_NAME}")
	return()
endif(NOT muParser_FOUND)

include_directories(${muParser_INCLUDE_DIRS})
link_directories(${muParser_LIBRARY_DIRS})

add_executable(cstool main.cc cross-section.cc
    compile-mat.cc plot-mat.cc ../common/colours.cc
	../common/parser.cc ../common/material.cc)

target_link_libraries(cstool common ${muParser_LIBRARIES})
