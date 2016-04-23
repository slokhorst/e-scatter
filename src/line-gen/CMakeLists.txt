project(e-scatter/line-gen)
include(FindPkgConfig)

find_package(Boost COMPONENTS program_options)
pkg_search_module(FFTW fftw3)
if(NOT (Boost_FOUND AND FFTW_FOUND))
	message(WARNING "missing libraries, not building ${PROJECT_NAME}")
	return()
endif(NOT (Boost_FOUND AND FFTW_FOUND))

include_directories(${Boost_INCLUDE_DIRS} ${FFTW_INCLUDE_DIRS})
add_executable(line-gen main.cc vertex_grid.cc ../cdsem/point3.cc ../cdsem/triangle.cc)
target_link_libraries(line-gen ${Boost_LIBRARIES} ${FFTW_LIBRARIES})
