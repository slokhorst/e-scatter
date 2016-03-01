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
 * @file src/cdsem/vector3.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SCATTER__CDSEM__VECTOR3__HEADER_INCLUDED
#define SCATTER__CDSEM__VECTOR3__HEADER_INCLUDED

inline __host__ __device__ float3 operator+(float3 a, float3 b);
inline __host__ __device__ float3 operator-(float3 a, float3 b);
inline __host__ __device__ float3 operator*(float3 a, float b);
inline __host__ __device__ float3 operator*(float a, float3 b);
inline __host__ __device__ float3 operator-(float3 a);
inline __host__ __device__ float3 cross_product(float3, float3);
inline __host__ __device__ float dot_product(float3, float3);
inline __host__ __device__ float norm2(float3);
inline __host__ __device__ float3 normalize(float3);

#include "vector3.cuinc"

#endif