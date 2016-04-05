project(e-scatter/sample-viewer)
include(FindPkgConfig)
pkg_search_module(SDL2 REQUIRED sdl2)
find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)
include_directories(${GLEW_INCLUDE_DIRS})
include_directories(${SDL2_INCLUDE_DIRS})
include_directories(${OPENGL_INCLUDE_DIR})
add_executable(sample-viewer main.cc sample_viewer.cc
	../common/archive.cc
	../cpl/triangle.cc ../cpl/text.cc ../cpl/vector3.cc ../cpl/parser.cc)
target_link_libraries(sample-viewer ${GLEW_LIBRARIES})
target_link_libraries(sample-viewer ${SDL2_LIBRARIES})
target_link_libraries(sample-viewer ${OPENGL_LIBRARIES})
