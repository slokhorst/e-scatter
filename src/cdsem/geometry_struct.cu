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
 * @file src/cdsem/geometry_struct.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "geometry_struct.cuh"
#include <algorithm>
#include <cfloat>
#include <functional>
#include <map>
#include <vector>
#include "cuda_mem_scope.cuh"
#include "cuda_safe_call.cuh"
#include "tri_box_overlap.h"

__host__ geometry_struct geometry_struct::malloc(float3 min, float3 max, float3 delta) {
	geometry_struct gstruct;
	gstruct.min = make_float3(std::min(min.x, max.x), std::min(min.y, max.y), std::min(min.z, max.z));
	gstruct.max = make_float3(std::max(min.x, max.x), std::max(min.y, max.y), std::max(min.z, max.z));
	gstruct.delta = delta;
	gstruct.dim.x = std::abs((max.x-min.x)/delta.x)+1;
	gstruct.dim.y = std::abs((max.y-min.y)/delta.y)+1;
	gstruct.dim.z = std::abs((max.z-min.z)/delta.z)+1;
	gstruct.dim.w = 0;
	cuda_safe_call(__FILE__, __LINE__, [&] {
		cudaMalloc(&gstruct.map_dev_p, gstruct.dim.x*gstruct.dim.y*gstruct.dim.z*sizeof(int));
		cuda_mem_scope<int>(gstruct.map_dev_p, gstruct.dim.x*gstruct.dim.y*gstruct.dim.z, [&](int* map_p) {
			for(int i = 0; i < gstruct.dim.x*gstruct.dim.y*gstruct.dim.z; i++)
				map_p[i] = -1;
		});
	});
	gstruct.in_dev_p = gstruct.out_dev_p = nullptr;
	gstruct.Ax_dev_p = gstruct.Ay_dev_p = gstruct.Az_dev_p = nullptr;
	gstruct.Bx_dev_p = gstruct.By_dev_p = gstruct.Bz_dev_p = nullptr;
	gstruct.Cx_dev_p = gstruct.Cy_dev_p = gstruct.Cz_dev_p = nullptr;
	gstruct.pitch = 0;
	return gstruct;
}

__host__ void geometry_struct::free(geometry_struct& gstruct) {
	cuda_safe_call(__FILE__, __LINE__, [&] {
		cudaFree(gstruct.map_dev_p);
	});
	gstruct.map_dev_p = nullptr;
	gstruct.clear();
}

__host__ void geometry_struct::clear() {
	dim.w = 0;
	cuda_safe_call(__FILE__, __LINE__, [&] {
		for(char* p : {in_dev_p, out_dev_p})
			cudaFree(p);
		for(float* p : {Ax_dev_p, Ay_dev_p, Az_dev_p})
			cudaFree(p);
		for(float* p : {Bx_dev_p, By_dev_p, Bz_dev_p})
			cudaFree(p);
		for(float* p : {Cx_dev_p, Cy_dev_p, Cz_dev_p})
			cudaFree(p);
		cuda_mem_scope<int>(map_dev_p, dim.x*dim.y*dim.z, [&](int* map_p) {
			for(int i = 0; i < dim.x*dim.y*dim.z; i++)
				map_p[i] = -1;
		});
	});
	in_dev_p = out_dev_p = nullptr;
	Ax_dev_p = Ay_dev_p = Az_dev_p = nullptr;
	Bx_dev_p = By_dev_p = Bz_dev_p = nullptr;
	Cx_dev_p = Cy_dev_p = Cz_dev_p = nullptr;
	pitch = 0;
}

__host__ int geometry_struct::push(const float3* A, const float3* B, const float3* C, int n, char in, char out) {
	if(n <= 0)
		return 0;
	if((in == NOP) && (out == NOP))
		return 0;
	if((in == VACUUM) && (out == VACUUM))
		return 0;
	if((in >= 0) && (out >= 0) && (in == out))
		return 0;
	struct triangle {
		char in, out;
		float3 A, B, C;
	};
	std::map<int,std::map<int,std::map<int,std::vector<triangle>>>> grid_map;
	cuda_safe_call(__FILE__, __LINE__, [&] {
		cuda_mem_scope<int>(map_dev_p, dim.x*dim.y*dim.z, [&](int* map_p) {
			std::map<int,int> index_map;
			for(int i = 0; i < dim.x*dim.y*dim.z; i++)
				if(map_p[i] >= 0)
					index_map[i] = map_p[i];
			if(index_map.empty())
				return;
			const int2 size =
				make_int2(dim.w,1+std::max_element(index_map.cbegin(), index_map.cend(),index_map.value_comp())->second);
			cuda_mem_scope<char>(in_dev_p, pitch, size, [&](char** in_p) {
			cuda_mem_scope<char>(out_dev_p, pitch, size, [&](char** out_p) {
				cuda_mem_scope<float>(Ax_dev_p, pitch, size, [&](float** Ax_p) {
				cuda_mem_scope<float>(Ay_dev_p, pitch, size, [&](float** Ay_p) {
				cuda_mem_scope<float>(Az_dev_p, pitch, size, [&](float** Az_p) {
					cuda_mem_scope<float>(Bx_dev_p, pitch, size, [&](float** Bx_p) {
					cuda_mem_scope<float>(By_dev_p, pitch, size, [&](float** By_p) {
					cuda_mem_scope<float>(Bz_dev_p, pitch, size, [&](float** Bz_p) {
						cuda_mem_scope<float>(Cx_dev_p, pitch, size, [&](float** Cx_p) {
						cuda_mem_scope<float>(Cy_dev_p, pitch, size, [&](float** Cy_p) {
						cuda_mem_scope<float>(Cz_dev_p, pitch, size, [&](float** Cz_p) {
							for(auto cit = index_map.cbegin(); cit != index_map.cend(); cit++) {
								const int gid = cit->first;
								const int x = gid%dim.x;
								const int y = (gid/dim.x)%dim.y;
								const int z = (gid/(dim.x*dim.y))%dim.z;
								const int j = cit->second;
								for(int i = 0; i < dim.w; i++)
									if((in_p[j][i] != NOP) || (out_p[j][i] != NOP)) {
										const float3 A = make_float3(Ax_p[j][i], Ay_p[j][i], Az_p[j][i]);
										const float3 B = make_float3(Bx_p[j][i], By_p[j][i], Bz_p[j][i]);
										const float3 C = make_float3(Cx_p[j][i], Cy_p[j][i], Bz_p[j][i]);
										grid_map[z][y][x].push_back({in_p[j][i], out_p[j][i], A, B, C});
									}
							}
						});
						});
						});
					});
					});
					});
				});
				});
				});
			});
			});
		});
	});
	int nt = 0;
	for(int i = 0; i < n; i++) {
		const int x1 = (std::min({A[i].x, B[i].x, C[i].x})-min.x)/delta.x;
		const int y1 = (std::min({A[i].y, B[i].y, C[i].y})-min.y)/delta.y;
		const int z1 = (std::min({A[i].z, B[i].z, C[i].z})-min.z)/delta.z;
		if((x1 > dim.x-1) || (y1 > dim.y-1) || (z1 > dim.z-1))
			continue;
		const int x2 = (std::max({A[i].x, B[i].x, C[i].x})-min.x)/delta.x;
		const int y2 = (std::max({A[i].y, B[i].y, C[i].y})-min.y)/delta.y;
		const int z2 = (std::max({A[i].z, B[i].z, C[i].z})-min.z)/delta.z;
		if((x2 < 0) || (y2 < 0) || (z2 < 0))
			continue;
		for(int z = std::max(0, z1-1); z <= std::min(dim.z-1, z2+1); z++)
		for(int y = std::max(0, y1-1); y <= std::min(dim.y-1, y2+1); y++)
		for(int x = std::max(0, x1-1); x <= std::min(dim.x-1, x2+1); x++) {
			double boxcenter[3];
			double boxhalfsize[3];
			double triverts[3][3];
			boxcenter[0] = min.x+(0.5+x)*delta.x;
			boxcenter[1] = min.y+(0.5+y)*delta.y;
			boxcenter[2] = min.z+(0.5+z)*delta.z;
			boxhalfsize[0] = (1.0+FLT_EPSILON)*delta.x/2;
			boxhalfsize[1] = (1.0+FLT_EPSILON)*delta.y/2;
			boxhalfsize[2] = (1.0+FLT_EPSILON)*delta.z/2;
			triverts[0][0] = A[i].x; triverts[0][1] = A[i].y; triverts[0][2] = A[i].z;
			triverts[1][0] = B[i].x; triverts[1][1] = B[i].y; triverts[1][2] = B[i].z;
			triverts[2][0] = C[i].x; triverts[2][1] = C[i].y; triverts[2][2] = C[i].z;
			if(triBoxOverlap(boxcenter, boxhalfsize, triverts) > 0) {
				grid_map[z][y][x].push_back({in, out, A[i], B[i], C[i]});
				nt++;
			}
		}
	}
	clear();
	int2 size = make_int2(0,0);
	for(const auto& z_slice : grid_map)
	for(const auto& y_slice : z_slice.second)
	for(const auto& x_slice : y_slice.second)
		if(x_slice.second.size() > 0) {
			if(x_slice.second.size() > static_cast<size_t>(dim.w))
				dim.w = x_slice.second.size();
			size.y++;
		}
	size.x = dim.w;
	if(size.y > 0) {
		const int _sizeof = std::max(sizeof(char), sizeof(float));
		cuda_safe_call(__FILE__, __LINE__, [&] {
			size_t _pitch;
			for(char** p : {&in_dev_p, &out_dev_p})
				cudaMallocPitch(p, &_pitch, dim.w*_sizeof, size.y);
			for(float** p : {&Ax_dev_p, &Ay_dev_p, &Az_dev_p})
				cudaMallocPitch(p, &_pitch, dim.w*_sizeof, size.y);
			for(float** p : {&Bx_dev_p, &By_dev_p, &Bz_dev_p})
				cudaMallocPitch(p, &_pitch, dim.w*_sizeof, size.y);
			for(float** p : {&Cx_dev_p, &Cy_dev_p, &Cz_dev_p})
				cudaMallocPitch(p, &_pitch, dim.w*_sizeof, size.y);
			pitch = _pitch;
		});
		cuda_safe_call(__FILE__, __LINE__, [&] {
			cuda_mem_scope<char>(in_dev_p, pitch*size.y, [&](char* in_p) {
				for(int i = 0; i < pitch*size.y; i++)
					in_p[i] = NOP;
			});
			cuda_mem_scope<char>(out_dev_p, pitch*size.y, [&](char* out_p) {
				for(int i = 0; i < pitch*size.y; i++)
					out_p[i] = NOP;
			});
			cuda_mem_scope<int>(map_dev_p, dim.x*dim.y*dim.z, [&](int* map_p) {
				for(int i = 0; i < dim.x*dim.y*dim.z; i++)
					map_p[i] = -1;
				cuda_mem_scope<char>(in_dev_p, pitch, size, [&](char** in_p) {
				cuda_mem_scope<char>(out_dev_p, pitch, size, [&](char** out_p) {
					cuda_mem_scope<float>(Ax_dev_p, pitch, size, [&](float** Ax_p) {
					cuda_mem_scope<float>(Ay_dev_p, pitch, size, [&](float** Ay_p) {
					cuda_mem_scope<float>(Az_dev_p, pitch, size, [&](float** Az_p) {
						cuda_mem_scope<float>(Bx_dev_p, pitch, size, [&](float** Bx_p) {
						cuda_mem_scope<float>(By_dev_p, pitch, size, [&](float** By_p) {
						cuda_mem_scope<float>(Bz_dev_p, pitch, size, [&](float** Bz_p) {
							cuda_mem_scope<float>(Cx_dev_p, pitch, size, [&](float** Cx_p) {
							cuda_mem_scope<float>(Cy_dev_p, pitch, size, [&](float** Cy_p) {
							cuda_mem_scope<float>(Cz_dev_p, pitch, size, [&](float** Cz_p) {
								int j = -1;
								for(const auto& z_slice : grid_map)
								for(const auto& y_slice : z_slice.second)
								for(const auto& x_slice : y_slice.second)
									if(!x_slice.second.empty()) {
										const int x = x_slice.first;
										const int y = y_slice.first;
										const int z = z_slice.first;
										map_p[x+y*dim.x+z*dim.y*dim.x] = ++j;
										for(size_t i = 0; i < x_slice.second.size(); i++) {
											in_p[j][i] = x_slice.second[i].in;
											out_p[j][i] = x_slice.second[i].out;
											Ax_p[j][i] = x_slice.second[i].A.x;
											Ay_p[j][i] = x_slice.second[i].A.y;
											Az_p[j][i] = x_slice.second[i].A.z;
											Bx_p[j][i] = x_slice.second[i].B.x;
											By_p[j][i] = x_slice.second[i].B.y;
											Bz_p[j][i] = x_slice.second[i].B.z;
											Cx_p[j][i] = x_slice.second[i].C.x;
											Cy_p[j][i] = x_slice.second[i].C.y;
											Cz_p[j][i] = x_slice.second[i].C.z;
										}
									}
							});
							});
							});
						});
						});
						});
					});
					});
					});
				});
				});
			});
		});
	}
	return nt;
}

__host__ int geometry_struct::count() const {
	int j = -1;
	cuda_safe_call(__FILE__, __LINE__, [&] {
		cuda_mem_scope<int>(map_dev_p, dim.x*dim.y*dim.z, [&](int* map_p) {
			for(int i = 0; i < dim.x*dim.y*dim.z; i++)
				j = std::max(map_p[i],j);
		});
	});
	return j+1;
}