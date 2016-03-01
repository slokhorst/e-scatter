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
 * @file src/cdsem/kernels.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SCATTER__CDSEM__KERNELS__HEADER_INCLUDED
#define SCATTER__CDSEM__KERNELS__HEADER_INCLUDED

#include <curand_kernel.h>
#include "geometry_struct.cuh"
#include "material_struct.cuh"
#include "particle_struct.cuh"

__global__ void __init_rand_state(curandState* rand_state_p, unsigned long long seed, int n);
__global__ void __init_trajectory(particle_struct pstruct, geometry_struct gstruct, material_struct mstruct);
__global__ void __update_trajectory(particle_struct pstruct, geometry_struct gstruct);
__global__ void __probe_isec_event(particle_struct pstruct, geometry_struct gstruct);
__global__ void __apply_isec_event(particle_struct pstruct, geometry_struct gstruct, material_struct mstruct, curandState* rand_state_dev_p);
__global__ void __probe_scatter_event(particle_struct pstruct, material_struct mstruct, curandState* rand_state_dev_p);
__global__ void __apply_elastic_event(particle_struct pstruct, material_struct mstruct, curandState* rand_state_dev_p);
__global__ void __apply_inelastic_event(particle_struct pstruct, material_struct mstruct, curandState* rand_state_dev_p);

#endif