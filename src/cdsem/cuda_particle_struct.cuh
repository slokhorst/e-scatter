/**
 * @file src/cdsem/cuda_particle_struct.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__CUDA_PARTICLE_STRUCT__HEADER_INCLUDED
#define eSCATTER__CDSEM__CUDA_PARTICLE_STRUCT__HEADER_INCLUDED

#include <functional>

struct cuda_particle_struct {
    __host__ static cuda_particle_struct create(int capacity);
    __host__ static void release(cuda_particle_struct&);

    __host__ void clear();
    __host__ void flush();
    __host__ int push(const float3* pos, const float3* dir, const float* K, const int* tag, int n);
    __host__ void for_each(int status, std::function<void(float3 pos, float3 dir, float K, int tag)>) const;

    enum status_enum : int {
        PENDING = 0,
        INELASTIC_EVENT,
        NEW_SECONDARY,
        ELASTIC_EVENT,
        INTERSECT_EVENT,
        NO_EVENT,
        DETECTED,
        TERMINATED
    };

    int capacity;
    int* status_dev_p;
    int* particle_idx_dev_p;
    int* particle_tag_dev_p;
    int* material_idx_dev_p;
    int* triangle_idx_dev_p;
    float* K_energy_dev_p;
    float* distance_dev_p;
    float* pos_x_dev_p; float* pos_y_dev_p; float* pos_z_dev_p;
    float* dir_x_dev_p; float* dir_y_dev_p; float* dir_z_dev_p;

private:
    __host__ __device__ cuda_particle_struct() = default;
};

#endif