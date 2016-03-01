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
 * @file src/cdsem/material_struct.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SCATTER__CDSEM__MATERIAL_STRUCT__HEADER_INCLUDED
#define SCATTER__CDSEM__MATERIAL_STRUCT__HEADER_INCLUDED

#include <functional>

//class MaterialStruct
//{
//    public:
//        __host__ static std::unique_ptr<MaterialStruct> create(int n);
//        __host__ __device__ MaterialStruct(MaterialStruct const &) = delete;
//
//        ~MaterialStruct() {} // free everything
//
//    private:
//        __host__ __device__ MaterialStruct() = default;
//};
//

struct material_struct {
	__host__ static material_struct malloc(int n);
	__host__ static void free(material_struct&);
	__host__ __device__ material_struct(const material_struct&) = default;
	__host__ __device__ material_struct& operator=(const material_struct&) = default;
	__host__ void set_inverse_mfp(int mid, std::function<void(const float* K, float* elastic_p, float* inelastic_p, int dim)>);
	__host__ void set_inverse_cdf(int mid, std::function<void(const float** K, const float** P, float** elastic_p, float** inelastic_p, int2 dim)>);
	__host__ void set_ionization(int mid, std::function<void(const float** K, const float** P, float** binding_p, int2 dim)>);
	int n;
	const int Kn = 64;
	const int Pn = 64;
	const float K1 = 1;
	const float K2 = 20e3;
	float* fermi_dev_p;
	float* barrier_dev_p;
	float* bandgap_dev_p;
	float* elastic_dev_p;
	float* inelastic_dev_p;
	float* ionization_dev_p;
	int pitch;
private:
	__host__ __device__ material_struct() = default;
};

#endif
