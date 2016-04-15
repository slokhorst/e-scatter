project(e-scatter/compile-mat)
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})
add_executable(compile-mat main.cc parser.cc
	../cdsem/material.cc 
	../cstool/cross-section.cc)
target_link_libraries(compile-mat common muparser ${Boost_LIBRARIES})
