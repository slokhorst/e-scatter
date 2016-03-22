/**
 * @file src/cdsem/cuda_particle_struct.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "cuda_particle_struct.cuh"
#include <vector>
#include <common/cuda_mem_scope.cuh>
#include <common/cuda_safe_call.cuh>
#include "cuda_geometry_struct.cuh"

__host__ cuda_particle_struct cuda_particle_struct::create(int capacity) {
	cuda_particle_struct pstruct;
	pstruct.capacity = capacity;
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cudaMalloc(&pstruct.status_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.pid_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.gid_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.tid_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.mid_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.tag_dev_p, capacity*sizeof(int));
		cudaMalloc(&pstruct.K_dev_p, capacity*sizeof(float));
		cudaMalloc(&pstruct.ds_dev_p, capacity*sizeof(float));
		for(float** p : {&pstruct.rx_dev_p, &pstruct.ry_dev_p, &pstruct.rz_dev_p})
			cudaMalloc(p, capacity*sizeof(float));
		for(float** p : {&pstruct.dx_dev_p, &pstruct.dy_dev_p, &pstruct.dz_dev_p})
			cudaMalloc(p, capacity*sizeof(float));
	});
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<int>(pstruct.status_dev_p, capacity, [&](int* status_p) {
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

__host__ void cuda_particle_struct::release(cuda_particle_struct& pstruct) {
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

__host__ void cuda_particle_struct::clear() {
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<int>(status_dev_p, capacity, [&](int* status_p) {
			for(int i = 0; i < capacity; i++)
				status_p[i] = TERMINATED;
		});
	});
}

__host__ void cuda_particle_struct::flush() {
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<int>(status_dev_p, capacity, [&](int* status_p) {
			for(int i = 0; i < capacity; i++)
				if(status_p[i] == DETECTED)
					status_p[i] = TERMINATED;
		});
	});
}

__host__ void cuda_particle_struct::for_each(int status, std::function<void(float3 pos, float3 dir, float K, int tag)> callback) const {
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<int>(status_dev_p, capacity, [&](int* status_p) {
		cuda_mem_scope<float>(rx_dev_p, capacity, [&](float* rx_p) {
		cuda_mem_scope<float>(ry_dev_p, capacity, [&](float* ry_p) {
		cuda_mem_scope<float>(rz_dev_p, capacity, [&](float* rz_p) {
		cuda_mem_scope<float>(dx_dev_p, capacity, [&](float* dx_p) {
		cuda_mem_scope<float>(dy_dev_p, capacity, [&](float* dy_p) {
		cuda_mem_scope<float>(dz_dev_p, capacity, [&](float* dz_p) {
		cuda_mem_scope<float>(K_dev_p, capacity, [&](float* K_p) {
		cuda_mem_scope<int>(tag_dev_p, capacity, [&](int* tag_p) {
			for(int i = 0; i < capacity; i++)
				if(status_p[i] == status) {
					float3 pos = make_float3(rx_p[i], ry_p[i], rz_p[i]);
					float3 dir = make_float3(dx_p[i], dy_p[i], dz_p[i]);
					callback(pos, dir, K_p[i], tag_p[i]);
				}
		}); /* tag_dev_p */
		}); /* rz_dev_p */
		}); /* ry_dev_p */
		}); /* rx_dev_p */
		}); /* dz_dev_p */
		}); /* dy_dev_p */
		}); /* dx_dev_p */
		}); /* K_dev_p */
		}); /* status_dev_p */
	});
}

__host__ int cuda_particle_struct::push(const float3* pos, const float3* dir, const float* K, const int* tag, int n) {
	if((pos == nullptr) || (dir == nullptr) || (K == nullptr) || (n <= 0))
		return 0;
	std::vector<int> index_vec;
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<int>(status_dev_p, capacity, [&](int* status_p) {
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
		cuda_mem_scope<int>(mid_dev_p, capacity, [&index_vec](int* mid_p) {
			for(int i : index_vec)
				mid_p[i] = cuda_geometry_struct::VACUUM;
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