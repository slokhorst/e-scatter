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

struct scatter_options {
	// flags to enable or disable specific losses
	bool acoustic_phonon_loss_flag = true;
	bool optical_phonon_loss_flag = true;
	bool atomic_recoil_loss_flag = true;

	// flag to enable or disable creation of secondary electrons
	bool generate_secondary_flag = true;

	// use a random instantaneous momentum of the secondary prior to collision?
	bool instantaneous_momentum_flag = true;

	// flag to enable or disable momentum conservation
	// note: without momentum conservation -> forward scattering of the incident electron is assumed
	bool momentum_conservation_flag = true;

	// interface related flags
	bool quantum_transmission_flag = true;
	bool interface_refraction_flag = true;
	bool interface_absorption_flag = true; // for compliance with Kieft & Bosch
};

__global__ void cuda_init_rand_state(curandState* rand_state_p, unsigned long long seed, int n, scatter_options opt);
__global__ void cuda_init_trajectory(cuda_particle_struct, cuda_geometry_struct, cuda_material_struct, curandState*, scatter_options opt);
__global__ void cuda_update_trajectory(cuda_particle_struct, cuda_geometry_struct, cuda_material_struct, scatter_options opt);
__global__ void cuda_intersection_event(cuda_particle_struct, cuda_geometry_struct, cuda_material_struct, curandState*, scatter_options opt);
__global__ void cuda_elastic_event(cuda_particle_struct, cuda_material_struct, curandState*, scatter_options opt);
__global__ void cuda_inelastic_event(cuda_particle_struct, cuda_material_struct, curandState*, scatter_options opt);

#endif