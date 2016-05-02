/**
 * @file src/cdsem/cuda_particle_struct.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "cuda_particle_struct.cuh"
#include <vector>
#include <common/cuda_mem_scope.cuh>
#include <common/cuda_safe_call.cuh>
#include "triangle.hh"

__host__ cuda_particle_struct cuda_particle_struct::create(int capacity) {
    cuda_particle_struct pstruct;
    pstruct.capacity = capacity;
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&pstruct.status_dev_p, capacity*sizeof(uint8_t));
        cudaMalloc(&pstruct.particle_idx_dev_p, capacity*sizeof(int));
        cudaMalloc(&pstruct.particle_tag_dev_p, capacity*sizeof(int));
        cudaMalloc(&pstruct.material_idx_dev_p, capacity*sizeof(int));
        cudaMalloc(&pstruct.K_energy_dev_p, capacity*sizeof(float));
        cudaMalloc(&pstruct.triangle_idx_dev_p, capacity*sizeof(int));
        cudaMalloc(&pstruct.distance_dev_p, capacity*sizeof(float));
        for(float** p : {&pstruct.pos_x_dev_p, &pstruct.pos_y_dev_p, &pstruct.pos_z_dev_p})
            cudaMalloc(p, capacity*sizeof(float));
        for(float** p : {&pstruct.dir_x_dev_p, &pstruct.dir_y_dev_p, &pstruct.dir_z_dev_p})
            cudaMalloc(p, capacity*sizeof(float));
    });

    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cuda_mem_scope<uint8_t>(pstruct.status_dev_p, capacity, [&](uint8_t* status_p) {
            for(int i = 0; i < capacity; i++)
                status_p[i] = TERMINATED;
        });
        cuda_mem_scope<int>(pstruct.particle_idx_dev_p, capacity, [&](int* pid_p) {
            for(int i = 0; i < capacity; i++)
                pid_p[i] = i;
        });
    });

    return pstruct;
};

__host__ void cuda_particle_struct::release(cuda_particle_struct& pstruct) {
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaFree(pstruct.status_dev_p);
        cudaFree(pstruct.particle_idx_dev_p);
        cudaFree(pstruct.triangle_idx_dev_p);
        cudaFree(pstruct.material_idx_dev_p);
        cudaFree(pstruct.particle_tag_dev_p);
        cudaFree(pstruct.K_energy_dev_p);
        cudaFree(pstruct.distance_dev_p);
        for(float* p : {pstruct.pos_x_dev_p, pstruct.pos_y_dev_p, pstruct.pos_z_dev_p})
            cudaFree(p);
        for(float* p : {pstruct.dir_x_dev_p, pstruct.dir_y_dev_p, pstruct.dir_z_dev_p})
            cudaFree(p);
    });
}

__host__ void cuda_particle_struct::clear() {
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cuda_mem_scope<uint8_t>(status_dev_p, capacity, [&](uint8_t* status_p) {
            for(int i = 0; i < capacity; i++)
                status_p[i] = TERMINATED;
        });
    });
}

__host__ void cuda_particle_struct::flush() {
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cuda_mem_scope<uint8_t>(status_dev_p, capacity, [&](uint8_t* status_p) {
            for(int i = 0; i < capacity; i++)
                if(status_p[i] == DETECTED)
                    status_p[i] = TERMINATED;
        });
    });
}

__host__ void cuda_particle_struct::for_each(uint8_t status, std::function<void(float3 pos, float3 dir, float K, int tag)> callback) const {
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cuda_mem_scope<uint8_t>(status_dev_p, capacity, [&](uint8_t* status_p) {
        cuda_mem_scope<float>(pos_x_dev_p, capacity, [&](float* pos_x_p) {
        cuda_mem_scope<float>(pos_y_dev_p, capacity, [&](float* pos_y_p) {
        cuda_mem_scope<float>(pos_z_dev_p, capacity, [&](float* pos_z_p) {
        cuda_mem_scope<float>(dir_x_dev_p, capacity, [&](float* dir_x_p) {
        cuda_mem_scope<float>(dir_y_dev_p, capacity, [&](float* dir_y_p) {
        cuda_mem_scope<float>(dir_z_dev_p, capacity, [&](float* dir_z_p) {
        cuda_mem_scope<float>(K_energy_dev_p, capacity, [&](float* K_energy_p) {
        cuda_mem_scope<int>(particle_tag_dev_p, capacity, [&](int* particle_tag_p) {
            for(int i = 0; i < capacity; i++)
                if(status_p[i] == status) {
                    float3 pos = make_float3(pos_x_p[i], pos_y_p[i], pos_z_p[i]);
                    float3 dir = make_float3(dir_x_p[i], dir_y_p[i], dir_z_p[i]);
                    callback(pos, dir, K_energy_p[i], particle_tag_p[i]);
                }
        }); // particle_tag_dev_p
        }); // pos_z_dev_p
        }); // pos_y_dev_p
        }); // pos_x_dev_p
        }); // dir_z_dev_p
        }); // dir_y_dev_p
        }); // dir_x_dev_p
        }); // K_energy_dev_p
        }); // status_dev_p
    });
}

__host__ int cuda_particle_struct::push(const float3* pos, const float3* dir, const float* K, const int* tag, int n) {
    if((pos == nullptr) || (dir == nullptr) || (K == nullptr) || (n <= 0))
        return 0;

    std::vector<int> free_index_vec;
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cuda_mem_scope<uint8_t>(status_dev_p, capacity, [&](uint8_t* status_p) {
            for(int i = 0; i < capacity; i++)
                if(free_index_vec.size() < static_cast<size_t>(n))
                    if(status_p[i] == TERMINATED) {
                        status_p[i] = NO_EVENT;
                        free_index_vec.push_back(i);
                    }
        });
    });

    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cuda_mem_scope<int>(triangle_idx_dev_p, capacity, [&free_index_vec](int* triangle_idx_p) {
            for(int i : free_index_vec)
                triangle_idx_p[i] = -1;
        });
        cuda_mem_scope<int>(material_idx_dev_p, capacity, [&free_index_vec](int* material_idx_p) {
            for(int i : free_index_vec)
                material_idx_p[i] = triangle::VACUUM;
        });
    });

    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cuda_mem_scope<int>(particle_tag_dev_p, capacity, [&free_index_vec,&tag](int* particle_tag_p) {
            for(size_t i = 0; i < free_index_vec.size(); i++)
                particle_tag_p[free_index_vec[i]] = tag[i];
        });
        cuda_mem_scope<float>(K_energy_dev_p, capacity, [&](float* K_energy_p) {
            for(size_t i = 0; i < free_index_vec.size(); i++)
                K_energy_p[free_index_vec[i]] = K[i];
        });
        cuda_mem_scope<float>(pos_x_dev_p, capacity, [&](float* pos_x_p) {
        cuda_mem_scope<float>(pos_y_dev_p, capacity, [&](float* pos_y_p) {
        cuda_mem_scope<float>(pos_z_dev_p, capacity, [&](float* pos_z_p) {
            for(size_t i = 0; i < free_index_vec.size(); i++) {
                pos_x_p[free_index_vec[i]] = pos[i].x;
                pos_y_p[free_index_vec[i]] = pos[i].y;
                pos_z_p[free_index_vec[i]] = pos[i].z;
            }
        });
        });
        });
        cuda_mem_scope<float>(dir_x_dev_p, capacity, [&](float* dir_x_p) {
        cuda_mem_scope<float>(dir_y_dev_p, capacity, [&](float* dir_y_p) {
        cuda_mem_scope<float>(dir_z_dev_p, capacity, [&](float* dir_z_p) {
            for(size_t i = 0; i < free_index_vec.size(); i++) {
                dir_x_p[free_index_vec[i]] = dir[i].x;
                dir_y_p[free_index_vec[i]] = dir[i].y;
                dir_z_p[free_index_vec[i]] = dir[i].z;
            }
        });
        });
        });
    });

    return free_index_vec.size();
}