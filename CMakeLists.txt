cmake_minimum_required(VERSION 2.8)

option(WITH_UTILS "Build all utilities." OFF)
option(WITH_TEST "Build all tests." OFF)
option(WITH_DOCS "Build documentation." OFF)
set(PROJECT_NAME eScatter-CUDA)
project(${PROJECT_NAME})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -Wall -ggdb")

find_package(CUDA REQUIRED)
list(APPEND CUDA_NVCC_FLAGS "-std=c++11 -O3 -prec-sqrt=false -Xptxas=-v --compiler-options -Wall")
list(APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_20,code=sm_20")
list(APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_52,code=sm_52")
set(CUDA_PROPAGATE_HOST_FLAGS off)

if(WITH_DOCS)
	find_package(Doxygen REQUIRED)
	add_custom_target(doc ALL ${DOXYGEN_EXECUTABLE} Doxyfile COMMENT "Generating documentation" VERBATIM)
endif(WITH_DOCS)

if (WITH_TEST)
	enable_testing()
	find_package(GTest REQUIRED)
	add_subdirectory(test)
endif(WITH_TEST)


include_directories(src 3rdparty/cub 3rdparty/GSL)
add_subdirectory(src)

