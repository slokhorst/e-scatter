project(e-scatter/sample-viewer)
include(FindPkgConfig)

pkg_search_module(SDL2 sdl2)
find_package(GLEW)
find_package(OpenGL)
if(NOT (SDL2_FOUND AND GLEW_FOUND AND OpenGL_FOUND))
	message(WARNING "missing libraries, not building ${PROJECT_NAME}")
	return()
endif(NOT (SDL2_FOUND AND GLEW_FOUND AND OpenGL_FOUND))

include_directories(${GLEW_INCLUDE_DIRS} ${SDL2_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR})
link_directories(${SDL2_LIBRARY_DIRS})

add_executable(sample-viewer main.cc sample_viewer.cc ../cdsem/point3.cc ../cdsem/triangle.cc)
target_link_libraries(sample-viewer common ${GLEW_LIBRARIES} ${SDL2_LIBRARIES} ${OPENGL_LIBRARIES})
