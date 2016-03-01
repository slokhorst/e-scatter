/*
 * Copyright 2014-2016 Thomas Verduin <T.Verduin@tudelft.nl>
 *
 * This program is part of the electron-matter interaction program (SCATTER).
 *
 * SCATTER is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

/**
 * @file src/cuda_mem_scope.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SCATTER__CUDA_MEM_SCOPE__HEADER_INCLUDED
#define SCATTER__CUDA_MEM_SCOPE__HEADER_INCLUDED

#include <array>
#include <functional>

template<typename T>
__host__ void cuda_mem_scope(T* dev_p, int dim, std::function<void(T* host_p)>);
template<typename T>
__host__ void cuda_mem_scope(T* dev_p, int pitch, int2 dim, std::function<void(T** host_p)>);
template<typename T>
__host__ void cuda_mem_scope(T* dev_p, int pitch, int height, int3 dim, std::function<void(T*** host_p)>);

#include "cuda_mem_scope.cuinc"

#endif
