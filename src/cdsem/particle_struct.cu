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
 * @file src/cdsem/particle_struct.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "particle_struct.cuh"
#include <vector>
#include "cuda_mem_scope.cuh"
#include "cuda_safe_call.cuh"
#include "geometry_struct.cuh"

__host__ particle_struct particle_struct::malloc(int capacity) {
	particle_struct pstruct;
	pstruct.capacity = capacity;
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cudaMalloc(&pstruct.status_dev_p, capacity*sizeof(unsigned char));
		cudaMalloc(&pstruct.pid_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.gid_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.tid_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.mid_dev_p, capacity*sizeof(char));
		cudaMalloc(&pstruct.tag_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.K_dev_p, capacity*sizeof(float));
		cudaMalloc(&pstruct.ds_dev_p, capacity*sizeof(float));
		for(float** p : {&pstruct.rx_dev_p, &pstruct.ry_dev_p, &pstruct.rz_dev_p})
			cudaMalloc(p, capacity*sizeof(float));
		for(float** p : {&pstruct.dx_dev_p, &pstruct.dy_dev_p, &pstruct.dz_dev_p})
			cudaMalloc(p, capacity*sizeof(float));
	});
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<unsigned char>(pstruct.status_dev_p, capacity, [&](unsigned char* status_p) {
			for(int i = 0; i < capacity; i++)
				status_p[i] = TERMINATED;
		});
		cuda_mem_scope<int>(pstruct.pid_dev_p, capacity, [&](int* pid_p) {
			for(int i = 0; i < capacity; i++)
				pid_p[i] = i;
		});
	});
	return pstruct;
};

__host__ void particle_struct::free(particle_struct& pstruct) {
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cudaFree(pstruct.status_dev_p);
		cudaFree(pstruct.pid_dev_p);
		cudaFree(pstruct.gid_dev_p);
		cudaFree(pstruct.tid_dev_p);
		cudaFree(pstruct.mid_dev_p);
		cudaFree(pstruct.tag_dev_p);
		cudaFree(pstruct.K_dev_p);
		cudaFree(pstruct.ds_dev_p);
		for(float* p : {pstruct.rx_dev_p, pstruct.ry_dev_p, pstruct.rz_dev_p})
			cudaFree(p);
		for(float* p : {pstruct.dx_dev_p, pstruct.dy_dev_p, pstruct.dz_dev_p})
			cudaFree(p);
	});
}

__host__ void particle_struct::clear() {
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<unsigned char>(status_dev_p, capacity, [&](unsigned char* status_p) {
			for(int i = 0; i < capacity; i++)
				status_p[i] = TERMINATED;
		});
	});
}

__host__ int particle_struct::push(const float3* pos, const float3* dir, const float* K, const int* tag, int n) {
	if(n <= 0)
		return 0;
	std::vector<int> index_vec;
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<unsigned char>(status_dev_p, capacity, [&](unsigned char* status_p) {
			for(int i = 0; i < capacity; i++)
				if(index_vec.size() < static_cast<size_t>(n))
					if(status_p[i] == TERMINATED) {
						status_p[i] = NO_EVENT;
						index_vec.push_back(i);
					}
		});
		cuda_mem_scope<int>(gid_dev_p, capacity, [&index_vec](int* gid_p) {
			for(int i : index_vec)
				gid_p[i] = -1;
		});
		cuda_mem_scope<char>(mid_dev_p, capacity, [&index_vec](char* mid_p) {
			for(int i : index_vec)
				mid_p[i] = geometry_struct::VACUUM;
		});
		cuda_mem_scope<int>(tag_dev_p, capacity, [&index_vec,&tag](int* tag_p) {
			for(size_t i = 0; i < index_vec.size(); i++)
				tag_p[index_vec[i]] = tag[i];
		});
		cuda_mem_scope<float>(K_dev_p, capacity, [&](float* K_p) {
			for(size_t i = 0; i < index_vec.size(); i++)
				K_p[index_vec[i]] = K[i];
		});
		cuda_mem_scope<float>(rx_dev_p, capacity, [&](float* rx_p) {
		cuda_mem_scope<float>(ry_dev_p, capacity, [&](float* ry_p) {
		cuda_mem_scope<float>(rz_dev_p, capacity, [&](float* rz_p) {
			for(size_t i = 0; i < index_vec.size(); i++) {
				rx_p[index_vec[i]] = pos[i].x;
				ry_p[index_vec[i]] = pos[i].y;
				rz_p[index_vec[i]] = pos[i].z;
			}
		});
		});
		});
		cuda_mem_scope<float>(dx_dev_p, capacity, [&](float* dx_p) {
		cuda_mem_scope<float>(dy_dev_p, capacity, [&](float* dy_p) {
		cuda_mem_scope<float>(dz_dev_p, capacity, [&](float* dz_p) {
			for(size_t i = 0; i < index_vec.size(); i++) {
				dx_p[index_vec[i]] = dir[i].x;
				dy_p[index_vec[i]] = dir[i].y;
				dz_p[index_vec[i]] = dir[i].z;
			}
		});
		});
		});
	});
	return index_vec.size();
}