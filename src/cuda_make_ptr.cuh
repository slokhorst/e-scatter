/*
 * Copyright 2014-2016 Thomas Verduin <T.Verduin@tudelft.nl>
 *
 * This program is part of the electron-matter interaction program (SCATTER).
 *
 * CPL is free software; you can redistribute it and/or modify
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
 * @file src/cuda_make_ptr.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SCATTER__CUDA_MAKE_PTR__HEADER_INCLUDED
#define SCATTER__CUDA_MAKE_PTR__HEADER_INCLUDED

template<typename T>
__host__ __device__ const T* cuda_make_ptr(const T* ptr, int pitch, int iy);
template<typename T>
__host__ __device__ const T* cuda_make_ptr(const T* ptr, int pitch, int height, int iy, int iz);

template<typename T>
__host__ __device__ T* cuda_make_ptr(T* ptr, int pitch, int iy, int iz);
template<typename T>
__host__ __device__ T* cuda_make_ptr(T* ptr, int pitch, int height, int iy, int iz);

#include "cuda_make_ptr.cuinc"

#endif
