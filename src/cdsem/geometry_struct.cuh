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
 * @file src/cdsem/geometry_struct.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SCATTER__CDSEM__GEOMETRY_STRUCT__HEADER_INCLUDED
#define SCATTER__CDSEM__GEOMETRY_STRUCT__HEADER_INCLUDED

struct geometry_struct {
	__host__ __device__ geometry_struct(const geometry_struct&) = default;
	__host__ __device__ geometry_struct& operator=(const geometry_struct&) = default;
	__host__ static geometry_struct malloc(float3 min, float3 max, float3 delta);
	__host__ static void free(geometry_struct&);
	__host__ void clear();
	__host__ int push(const float3* A, const float3* B, const float3* C, int n, char in, char out);
	__host__ int count() const;
	enum id_enum : char {
		NOP = -128,
		TERMINATOR, // TODO
		DETECTOR,
		DETECTOR_LT50,
		DETECTOR_GE50,
		VACUUM,
		MIRROR, // TODO
	};
	float3 min, max, delta;
	int4 dim;
	int* map_dev_p;
	char* in_dev_p; char* out_dev_p;
	float* Ax_dev_p; float* Ay_dev_p; float* Az_dev_p;
	float* Bx_dev_p; float* By_dev_p; float* Bz_dev_p;
	float* Cx_dev_p; float* Cy_dev_p; float* Cz_dev_p;
	int pitch;
private:
	__host__ __device__ geometry_struct() = default;
};

#endif