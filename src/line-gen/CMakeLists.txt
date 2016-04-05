project(e-scatter/line-gen)
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})
add_executable(line-gen main.cc vertex_grid.cc
	../cpl/text.cc ../cpl/triangle.cc ../cpl/vector3.cc ../cpl/parser.cc)
target_link_libraries(line-gen ${Boost_LIBRARIES})
target_link_libraries(line-gen fftw3)
