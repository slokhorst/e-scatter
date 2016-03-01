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
 * @file src/cdsem/cdsem.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <cmath>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <vector>
#include <curand_kernel.h>
#include <cub/cub.cuh>
#include "cuda_mem_scope.cuh"
#include "cuda_safe_call.cuh"
#include "interpolate.h"
#include "kernels.cuh"
#include "geometry_struct.cuh"
#include "material_struct.cuh"
#include "particle_struct.cuh"

const int warp_size = 32;
const int block_size = warp_size*4;
int block_count;
int* radix_index_dev_p;
void* radix_dump_dev_p;
void* radix_temp_dev_p;
size_t radix_temp_size;
curandState* rand_state_dev_p;

void load_tri_file(geometry_struct& gstruct, const std::string& file) {
	std::ifstream ifs(file);
	struct triangle {
		float3 A, B, C;
	};
	std::map<std::pair<int,int>,std::vector<triangle>> tri_map;
	while(!ifs.eof()) {
		int in, out;
		ifs >> in >> out;
		float3 A, B, C;
		ifs >> A.x >> A.y >> A.z;
		ifs >> B.x >> B.y >> B.z;
		ifs >> C.x >> C.y >> C.z;
		if(ifs.eof())
			break;
		tri_map[std::make_pair(in, out)].push_back({A, B, C});
	}
	ifs.close();
	for(auto cit = tri_map.cbegin(); cit != tri_map.cend(); cit++) {
		float3* A = new float3[cit->second.size()];
		float3* B = new float3[cit->second.size()];
		float3* C = new float3[cit->second.size()];
		for(size_t i = 0; i < cit->second.size(); i++) {
			A[i] = cit->second[i].A;
			B[i] = cit->second[i].B;
			C[i] = cit->second[i].C;
		}
		std::clog << gstruct.push(A, B, C, cit->second.size(), cit->first.first, cit->first.second) << std::endl;
		delete[] A;
		delete[] B;
		delete[] C;
	}
}

void load_mat_file(int i, material_struct& mstruct, const std::string& file) {
	std::vector<std::vector<std::vector<double>>> data_cube;
	int n = 0;
	std::ifstream ifs(file);
	while(!ifs.eof()) {
		std::string line;
		std::getline(ifs, line);
		if(n%2 == 0)
			if(!line.empty()) {
				data_cube.push_back(std::vector<std::vector<double>>());
				n++;
			}
		if(n%2 == 1) {
			if(!line.empty()) {
				data_cube.back().push_back(std::vector<double>());
				std::stringstream ss;
				ss << line;
				while(!ss.eof()) {
					float num;
					ss >> num;
					data_cube.back().back().push_back(num);
				}
			}
			else
				n++;
		}
	}
	ifs.close();
	const double fermi_energy = data_cube.at(0).at(0).at(0);
	const double barrier = data_cube.at(0).at(0).at(1);
	const double bandgap = data_cube.at(0).at(0).at(2);
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<float>(mstruct.fermi_dev_p, mstruct.n, [&](float* fermi_p) {
		cuda_mem_scope<float>(mstruct.barrier_dev_p, mstruct.n, [&](float* barrier_p) {
		cuda_mem_scope<float>(mstruct.bandgap_dev_p, mstruct.n, [&](float* bandgap_p) {
			fermi_p[i] = fermi_energy;
			barrier_p[i] = barrier;
			bandgap_p[i] = bandgap;
		});
		});
		});
	});
	std::map<double,double> elastic_imfp_map;
	std::map<double,std::map<double,double>> elastic_icdf_map;
	for(size_t i = 1; i < data_cube.at(1).size(); i++) {
		const double log_K = std::log(data_cube.at(1).at(i).at(0));
		const double imfp = 1.0/data_cube.at(1).at(i).at(1);
		elastic_imfp_map[log_K] = std::log(imfp);
		for(size_t j = 0; j < data_cube.at(1).at(0).size(); j++) {
			const double P = data_cube.at(1).at(0).at(j);
			const double theta = data_cube.at(1).at(i).at(2+j);
			elastic_icdf_map[log_K][P] = theta;
		}
	}
	std::map<double,double> inelastic_imfp_map;
	std::map<double,std::map<double,double>> inelastic_icdf_map;
	for(size_t i = 1; i < data_cube.at(2).size(); i++) {
		const double log_K = std::log(data_cube.at(2).at(i).at(0));
		const double imfp = 1.0/data_cube.at(2).at(i).at(1);
		inelastic_imfp_map[log_K] = std::log(imfp);
		for(size_t j = 0; j < data_cube.at(2).at(0).size(); j++) {
			const double P = data_cube.at(2).at(0).at(j);
			const double omega0 = data_cube.at(2).at(i).at(2+j);
			inelastic_icdf_map[log_K][P] = omega0;
		}
	}
	std::map<double,std::map<double,double>> ionization_map;
	for(size_t i = 2; i < data_cube.at(3).size(); i++) {
		const double log_K = std::log(data_cube.at(3).at(i).at(0));
		for(size_t j = 0; j < data_cube.at(3).at(1).size(); j++) {
			const double P = data_cube.at(3).at(1).at(j);
			const double B = data_cube.at(3).at(i).at(1+j);
			ionization_map[log_K][P] = B;
		}
	}
	mstruct.set_inverse_mfp(i, [&](const float* K, float* elastic_p, float* inelastic_p, int dim) {
		for(int i = 0; i < dim; i++) {
			elastic_p[i] = std::exp(interpolate(elastic_imfp_map, std::log(K[i])));
			inelastic_p[i] = std::exp(interpolate(inelastic_imfp_map, std::log(K[i])));
		}
	});
	mstruct.set_inverse_cdf(i, [&](const float** K, const float** P, float** elastic_p, float** inelastic_p, int2 dim) {
		for(int y = 0; y < dim.y; y++)
		for(int x = 0; x < dim.x; x++) {
			elastic_p[y][x] = std::cos(interpolate(elastic_icdf_map, std::log(K[y][x]), P[y][x]));
			inelastic_p[y][x] = interpolate(inelastic_icdf_map, std::log(K[y][x]), P[y][x]);
		}
	});
	cuda_mem_scope<float>(cuda_make_ptr<float>(mstruct.ionization_dev_p, mstruct.pitch, mstruct.Pn+1, 0, i), mstruct.Kn, [&](float* B) {
		for(int i = 0; i < mstruct.Kn; i++)
			B[i] = -1;
		for(size_t i = 0; i < data_cube.at(3).at(0).size(); i++)
			B[i] = data_cube.at(3).at(0).at(i);
	});
	mstruct.set_ionization(i, [&](const float** K, const float** P, float** binding_p, int2 dim) {
		for(int y = 0; y < dim.y; y++)
		for(int x = 0; x < dim.x; x++) {
			const double log_K = std::log(K[y][x]);
			double B = 0;
			if(log_K >= ionization_map.cbegin()->first) {
				auto cit = std::prev(ionization_map.upper_bound(log_K));
				if(P[y][x] >= cit->second.cbegin()->first)
					B = std::prev(cit->second.upper_bound(P[y][x]))->second;
			}
			binding_p[y][x] = B;
		}
	});
}

int main(const int argc, char* argv[]) {
	cuda_safe_call(__FILE__, __LINE__, [&] {
		cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
	});

	const float3 min = make_float3(-751, -751, -1999);
	const float3 max = make_float3(751, 751, 35);
	const float3 delta = make_float3(5, 5, 5);
	geometry_struct gstruct = geometry_struct::malloc(min, max, delta);
	load_tri_file(gstruct, "dat/rough-line.tri");
	const float3 size = make_float3(max.x-min.x, max.y-min.y, max.z-min.z);
	{ /* add detector (top)*/
		std::vector<float3> A_vec, B_vec, C_vec;
		A_vec.push_back(make_float3(-size.x, -size.y, 33));
		B_vec.push_back(make_float3(size.x, -size.y, 33));
		C_vec.push_back(make_float3(size.x, size.y, 33));
		A_vec.push_back(make_float3(-size.x, -size.y, 33));
		B_vec.push_back(make_float3(size.x, size.y, 33));
		C_vec.push_back(make_float3(-size.x, size.y, 33));
		std::clog << gstruct.push(A_vec.data(), B_vec.data(), C_vec.data(), 2, geometry_struct::NOP, geometry_struct::DETECTOR_LT50) << std::endl;
	}
	{ /* add detector (left)*/
		std::vector<float3> A_vec, B_vec, C_vec;
		A_vec.push_back(make_float3(0.1-size.x/2, -size.y, 0));
		B_vec.push_back(make_float3(0.1-size.x/2, -size.y, 33));
		C_vec.push_back(make_float3(0.1-size.x/2, size.y, 33));
		A_vec.push_back(make_float3(0.1-size.x/2, -size.y, 0));
		B_vec.push_back(make_float3(0.1-size.x/2, size.y, 33));
		C_vec.push_back(make_float3(0.1-size.x/2, size.y, 0));
		std::clog << gstruct.push(A_vec.data(), B_vec.data(), C_vec.data(), 2, geometry_struct::NOP, geometry_struct::DETECTOR_LT50) << std::endl;
	}
	{ /* add detector (right)*/
		std::vector<float3> A_vec, B_vec, C_vec;
		A_vec.push_back(make_float3(-0.1+size.x/2, -size.y, 33));
		B_vec.push_back(make_float3(-0.1+size.x/2, -size.y, 0));
		C_vec.push_back(make_float3(-0.1+size.x/2, size.y, 33));
		A_vec.push_back(make_float3(-0.1+size.x/2, -size.y, 0));
		B_vec.push_back(make_float3(-0.1+size.x/2, size.y, 0));
		C_vec.push_back(make_float3(-0.1+size.x/2, size.y, 33));
		std::clog << gstruct.push(A_vec.data(), B_vec.data(), C_vec.data(), 2, geometry_struct::NOP, geometry_struct::DETECTOR_LT50) << std::endl;
	}
	{ /* add detector (front)*/
		std::vector<float3> A_vec, B_vec, C_vec;
		A_vec.push_back(make_float3(-size.x, 0.1-size.y/2, 33));
		B_vec.push_back(make_float3(-size.x, 0.1-size.y/2, 0));
		C_vec.push_back(make_float3(size.x, 0.1-size.y/2, 0));
		A_vec.push_back(make_float3(-size.x, 0.1-size.y/2, 33));
		B_vec.push_back(make_float3(size.x, 0.1-size.y/2, 0));
		C_vec.push_back(make_float3(size.x, 0.1-size.y/2, 33));
		std::clog << gstruct.push(A_vec.data(), B_vec.data(), C_vec.data(), 2, geometry_struct::NOP, geometry_struct::DETECTOR_LT50) << std::endl;
	}
	{ /* add detector (back)*/
		std::vector<float3> A_vec, B_vec, C_vec;
		A_vec.push_back(make_float3(-size.x, -0.1+size.y/2, 33));
		B_vec.push_back(make_float3(size.x, -0.1+size.y/2, 0));
		C_vec.push_back(make_float3(-size.x, -0.1+size.y/2, 0));
		A_vec.push_back(make_float3(-size.x, -0.1+size.y/2, 33));
		B_vec.push_back(make_float3(size.x, -0.1+size.y/2, 33));
		C_vec.push_back(make_float3(size.x, -0.1+size.y/2, 0));
		std::clog << gstruct.push(A_vec.data(), B_vec.data(), C_vec.data(), 2, geometry_struct::NOP, geometry_struct::DETECTOR_LT50) << std::endl;
	}

	material_struct mstruct = material_struct::malloc(2);
	load_mat_file(0, mstruct, "dat/silicon.mat");
	load_mat_file(1, mstruct, "dat/pmma.mat");

	particle_struct pstruct = particle_struct::malloc(6250*32);
	block_count = 1+pstruct.capacity/block_size;
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cudaMalloc(&radix_index_dev_p, pstruct.capacity*sizeof(int));
		cuda_mem_scope<int>(radix_index_dev_p, pstruct.capacity, [&](int* index_p) {
			for(int i = 0; i < pstruct.capacity; i++)
				index_p[i] = i;
		});
		cudaMalloc(&radix_dump_dev_p, pstruct.capacity*sizeof(unsigned char));
		cub::DeviceRadixSort::SortPairs<unsigned char,int>(
			nullptr, radix_temp_size,
			pstruct.status_dev_p, static_cast<unsigned char*>(radix_dump_dev_p),
			radix_index_dev_p, pstruct.pid_dev_p,
			pstruct.capacity
		);
		radix_temp_size++;
		cudaMalloc(&radix_temp_dev_p, radix_temp_size);
		cudaMalloc(&rand_state_dev_p, pstruct.capacity*sizeof(curandState));
		__init_rand_state<<<block_count,block_size>>>(
			rand_state_dev_p, time(nullptr), pstruct.capacity
		);
	});

	std::map<int,std::pair<int,int>> tag_map;

	std::ifstream ifs("dat/exposure.pri");
	std::vector<float3> pos_vec;
	std::vector<float3> dir_vec;
	std::vector<float> K_vec;
	std::vector<int> tag_vec;
	while(!ifs.eof()) {
		int tag = tag_map.size();
		int tag_x, tag_y;
		ifs >> tag_x >> tag_y;
		float K;
		ifs >> K;
		float3 pos;
		ifs >> pos.x >> pos.y >> pos.z;
		float3 dir;
		ifs >> dir.x >> dir.y >> dir.z;
		if(ifs.eof())
			break;
		pos_vec.push_back(pos);
		dir_vec.push_back(dir);
		K_vec.push_back(K);
		tag_vec.push_back(tag);
		tag_map[tag] = std::make_pair(tag_x, tag_y);
	}
	ifs.close();

	const int np = tag_map.size();
	const int frame_size = 31;
	int batch_size = 3500;

	std::map<int,int> detector_map;
	int batch_index = 0;
	int running_count = batch_size;
	const std::chrono::high_resolution_clock::time_point chrono_start = std::chrono::high_resolution_clock::now();
	while(running_count > 0) {
		batch_size = std::min(np-batch_index, batch_size);
		if(batch_size > 0)
			batch_index += pstruct.push(pos_vec.data()+batch_index, dir_vec.data()+batch_index, K_vec.data()+batch_index, tag_vec.data()+batch_index, batch_size);
		cuda_safe_call(__FILE__, __LINE__, [&]() {
			for(int i = 0; i < frame_size; i++) {
				__init_trajectory<<<block_count,block_size>>>(
					pstruct, gstruct, mstruct
				);
				__probe_isec_event<<<block_count,block_size>>>(
					pstruct, gstruct
				);
				__probe_scatter_event<<<block_count,block_size>>>(
					pstruct, mstruct, rand_state_dev_p
				);
				__update_trajectory<<<block_count,block_size>>>(
					pstruct, gstruct
				);
				cub::DeviceRadixSort::SortPairs<unsigned char,int>(
					radix_temp_dev_p, radix_temp_size,
					pstruct.status_dev_p, static_cast<unsigned char*>(radix_dump_dev_p),
					radix_index_dev_p, pstruct.pid_dev_p,
					pstruct.capacity, 0, 4
				);
				__apply_isec_event<<<block_count,block_size>>>(
					pstruct, gstruct, mstruct, rand_state_dev_p
				);
				__apply_elastic_event<<<block_count,block_size>>>(
					pstruct, mstruct, rand_state_dev_p
				);
				__apply_inelastic_event<<<block_count,block_size>>>(
					pstruct, mstruct, rand_state_dev_p
				);
				cudaDeviceSynchronize();
			}
		});
		running_count = 0;
		int terminated_count = 0;
		cuda_safe_call(__FILE__, __LINE__, [&]() {
			cuda_mem_scope<unsigned char>(pstruct.status_dev_p, pstruct.capacity, [&](unsigned char* status_p) {
			cuda_mem_scope<int>(pstruct.tag_dev_p, pstruct.capacity, [&](int* tag_p) {
				for(int i = 0; i < pstruct.capacity; i++)
					switch(status_p[i]) {
						case particle_struct::INELASTIC_EVENT:
						case particle_struct::PENDING:
						case particle_struct::ELASTIC_EVENT:
						case particle_struct::ISEC_EVENT:
						case particle_struct::GRID_X_EVENT:
						case particle_struct::GRID_Y_EVENT:
						case particle_struct::GRID_Z_EVENT:
						case particle_struct::NEW_SECONDARY:
							running_count++;
						break;
						case particle_struct::TERMINATED:
							terminated_count++;
						break;
						case particle_struct::DETECTED:
							terminated_count++;
							status_p[i] = particle_struct::TERMINATED;
							detector_map[tag_p[i]]++;
						break;
						default:
						break;
					}
			}); // tag_dev_p
			}); // status_dev_p
		});
		std::clog << "\r" << 100*batch_index/(np-1) << "% ";
		std::clog << running_count << " ";
		std::clog << terminated_count << " ";
		std::clog << std::flush;
	}
	const std::chrono::high_resolution_clock::time_point chrono_stop = std::chrono::high_resolution_clock::now();
	std::clog << std::endl;
	std::clog << "time lapse=" << std::chrono::duration_cast<std::chrono::milliseconds>(chrono_stop-chrono_start).count() << "ms" << std::endl;

	for(auto cit = detector_map.cbegin(); cit != detector_map.cend(); cit++) {
		std::cout << tag_map[cit->first].first << " ";
		std::cout << tag_map[cit->first].second << " ";
		std::cout << cit->second << std::endl;
	}

	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cudaFree(radix_index_dev_p);
		cudaFree(radix_dump_dev_p);
		cudaFree(radix_temp_dev_p);
		cudaFree(rand_state_dev_p);
	});
	geometry_struct::free(gstruct);
	material_struct::free(mstruct);
	particle_struct::free(pstruct);

	return EXIT_SUCCESS;
}