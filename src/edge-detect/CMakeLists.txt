project(e-scatter/edge-detect)

find_package(Boost COMPONENTS program_options iostreams)
find_package(ImageMagick COMPONENTS Magick++)
find_package(OpenMP)
if(NOT (Boost_FOUND AND ImageMagick_FOUND AND OpenMP_FOUND))
	message(WARNING "missing libraries, not building ${PROJECT_NAME}")
	return()
endif(NOT (Boost_FOUND AND ImageMagick_FOUND AND OpenMP_FOUND))

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
add_definitions(-DMAGICKCORE_QUANTUM_DEPTH=16)
add_definitions(-DMAGICKCORE_HDRI_ENABLE=0)

include_directories(${Boost_INCLUDE_DIRS} ${ImageMagick_INCLUDE_DIRS})
add_executable(edge-detect main.cc float_image.cc gauss2_model.cc matrix.cc)
target_link_libraries(edge-detect ${Boost_LIBRARIES} ${ImageMagick_LIBRARIES})
