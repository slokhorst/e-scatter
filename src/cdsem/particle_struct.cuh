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
 * @file src/cdsem/particle_struct.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SCATTER__CDSEM__PARTICLE_STRUCT__HEADER_INCLUDED
#define SCATTER__CDSEM__PARTICLE_STRUCT__HEADER_INCLUDED

struct particle_struct {
	__host__ __device__ particle_struct(const particle_struct&) = default;
	__host__ __device__ particle_struct& operator=(const particle_struct&) = default;
	__host__ static particle_struct malloc(int capacity);
	__host__ static void free(particle_struct&);
	__host__ void clear();
	__host__ int push(const float3* pos, const float3* dir, const float* K, const int* tag, int n);
	enum status_enum : char {
		PENDING = 0,
		INELASTIC_EVENT,
		NEW_SECONDARY,
		ELASTIC_EVENT,
		ISEC_EVENT,
		GRID_X_EVENT,
		GRID_Y_EVENT,
		GRID_Z_EVENT,
		NO_EVENT,
		DETECTED,
		TERMINATED
	};
	int capacity;
	unsigned char* status_dev_p;
	int* pid_dev_p;
	int* gid_dev_p;
	int* tid_dev_p;
	char* mid_dev_p;
	int* tag_dev_p;  // 
	float* K_dev_p;  // kyn energy
	float* ds_dev_p; // attenuation length
	float* rx_dev_p; float* ry_dev_p; float* rz_dev_p;
	float* dx_dev_p; float* dy_dev_p; float* dz_dev_p;
private:
	__host__ __device__ particle_struct() = default;
};

#endif
