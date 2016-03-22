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
	__host__ void for_each(int status, std::function<void(float3 pos, float3 dir, float K, int tag)>) const;
	__host__ int push(const float3* pos, const float3* dir, const float* K, const int* tag, int n);
	enum status_enum : int {
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
	int* status_dev_p;
	int* pid_dev_p;
	int* gid_dev_p;
	int* tid_dev_p;
	int* mid_dev_p;
	int* tag_dev_p;
	float* K_dev_p;
	float* ds_dev_p;
	float* rx_dev_p; float* ry_dev_p; float* rz_dev_p;
	float* dx_dev_p; float* dy_dev_p; float* dz_dev_p;
private:
	__host__ __device__ cuda_particle_struct() = default;
};

#endif