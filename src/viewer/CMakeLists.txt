project(e-scatter/viewer)
include(FindPkgConfig)

pkg_search_module(SDL2 sdl2)
if(NOT SDL2_FOUND)
	message(WARNING "missing libraries, not building ${PROJECT_NAME}")
	return()
endif(NOT SDL2_FOUND)

include_directories(${SDL2_INCLUDE_DIRS})
link_directories(${SDL2_LIBRARY_DIRS})
add_executable(viewer main.cc)
target_link_libraries(viewer common ${SDL2_LIBRARIES})
