#  Copyright 2014-2016 Thomas Verduin <T.Verduin@tudelft.nl>
#
#  This file is part of the electron-matter interaction program (SCATTER).
#
#  CPL is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.

project(scatter)
cmake_minimum_required(VERSION 2.8)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -Wall -ggdb")
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} --compiler-options ${OpenMP_CXX_FLAGS}")
find_package(CUDA REQUIRED)
list(APPEND CUDA_NVCC_FLAGS "-std=c++11 -O3 -Xptxas=-v --compiler-options -Wall")
list(APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_20,code=sm_20")
list(APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_52,code=sm_52")
set(CUDA_PROPAGATE_HOST_FLAGS off)
find_package(Doxygen)
if(DOXYGEN_FOUND)
    add_custom_target(doc ALL ${DOXYGEN_EXECUTABLE} Doxyfile COMMENT "Generating documentation" VERBATIM)
endif(DOXYGEN_FOUND)
include_directories(3rdparty/cub-1.4.1)
include_directories(src)
add_subdirectory(src)
