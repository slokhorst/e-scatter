/**
 * @file src/cdsem/cuda_kernels.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__CUDA_KERNELS__HEADER_INCLUDED
#define eSCATTER__CDSEM__CUDA_KERNELS__HEADER_INCLUDED

#include <curand_kernel.h>
#include "cuda_geometry_struct.cuh"
#include "cuda_material_struct.cuh"
#include "cuda_particle_struct.cuh"

__global__ void cuda_init_rand_state(curandState* rand_state_p, unsigned long long seed, int n);
__global__ void cuda_init_trajectory(cuda_particle_struct, cuda_geometry_struct, cuda_material_struct, curandState*);
__global__ void cuda_update_trajectory(cuda_particle_struct, cuda_geometry_struct, cuda_material_struct);
__global__ void cuda_apply_intersection_event(cuda_particle_struct, cuda_geometry_struct, cuda_material_struct, curandState*);
__global__ void cuda_apply_elastic_event(cuda_particle_struct, cuda_material_struct, curandState*);
__global__ void cuda_apply_inelastic_event(cuda_particle_struct, cuda_material_struct, curandState*);

#endif