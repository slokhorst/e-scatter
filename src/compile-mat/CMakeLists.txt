project(e-scatter/compile-mat)
include(FindPkgConfig)

find_package(Boost COMPONENTS program_options)
pkg_search_module(muParser muparser)
if(NOT (Boost_FOUND AND muParser_FOUND))
	message(WARNING "missing libraries, not building ${PROJECT_NAME}")
	return()
endif(NOT (Boost_FOUND AND muParser_FOUND))

include_directories(${Boost_INCLUDE_DIRS} ${muParser_INCLUDE_DIRS})
add_executable(compile-mat main.cc parser.cc
	../cdsem/material.cc 
	../cstool/cross-section.cc)
target_link_libraries(compile-mat common ${muParser_LIBRARIES} ${Boost_LIBRARIES})
