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
 * @file src/cdsem/material_struct.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "material_struct.cuh"
#include "cuda_make_ptr.cuh"
#include "cuda_mem_scope.cuh"
#include "cuda_safe_call.cuh"

__host__ material_struct material_struct::malloc(int n) {
	material_struct mstruct;
	mstruct.n = n;
	cuda_safe_call(__FILE__, __LINE__, [&] {
		size_t _pitch;
		cudaMalloc(&mstruct.fermi_dev_p, n*sizeof(float));
		cudaMalloc(&mstruct.barrier_dev_p, n*sizeof(float));
		cudaMalloc(&mstruct.bandgap_dev_p, n*sizeof(float));
		cudaMallocPitch(&mstruct.elastic_dev_p, &_pitch, mstruct.Kn*sizeof(float), (mstruct.Pn+1)*n);
		cudaMallocPitch(&mstruct.inelastic_dev_p, &_pitch, mstruct.Kn*sizeof(float), (mstruct.Pn+1)*n);
		cudaMallocPitch(&mstruct.ionization_dev_p, &_pitch, mstruct.Kn*sizeof(float), (mstruct.Pn+1)*n);
		mstruct.pitch = _pitch;
	});
	return mstruct;
};

__host__ void material_struct::free(material_struct& mstruct) {
	cuda_safe_call(__FILE__, __LINE__, [&] {
		cudaFree(mstruct.fermi_dev_p);
		cudaFree(mstruct.barrier_dev_p);
		cudaFree(mstruct.bandgap_dev_p);
		cudaFree(mstruct.elastic_dev_p);
		cudaFree(mstruct.inelastic_dev_p);
		cudaFree(mstruct.ionization_dev_p);
	});
}

__host__ void material_struct::set_inverse_mfp(int mid, std::function<void(const float* K, float* elastic_p, float* inelastic_p, int dim)> callback) {
	float* elastic_imfp_dev_p = cuda_make_ptr<float>(elastic_dev_p, pitch, Pn+1, 0, mid);
	float* inelastic_imfp_dev_p = cuda_make_ptr<float>(inelastic_dev_p, pitch, Pn+1, 0, mid);
	cuda_safe_call(__FILE__, __LINE__, [&] {
		cuda_mem_scope<float>(elastic_imfp_dev_p, Kn, [&](float* elastic_p) {
		cuda_mem_scope<float>(inelastic_imfp_dev_p, Kn, [&](float* inelastic_p) {
			float* K = new float[Kn];
			for(int i = 0; i < Kn; i++) {
				K[i] = std::exp(std::log(K1)+std::log(K2/K1)*static_cast<double>(i)/(Kn-1));
				elastic_p[i] = std::exp(elastic_p[i]);
				inelastic_p[i] = std::exp(inelastic_p[i]);
			}
			callback(K, elastic_p, inelastic_p, Kn);
			for(int i = 0; i < Kn; i++) {
				elastic_p[i] = std::log(elastic_p[i]);
				inelastic_p[i] = std::log(inelastic_p[i]);
			}
			delete[] K;
		});
		});
	});
}

__host__ void material_struct::set_inverse_cdf(int mid, std::function<void(const float** K, const float** P, float** elastic_p, float** inelastic_p, int2 dim)> callback) {
	float* elastic_icdf_dev_p = cuda_make_ptr(elastic_dev_p, pitch, Pn+1, 1, mid);
	float* inelastic_icdf_dev_p = cuda_make_ptr(inelastic_dev_p, pitch, Pn+1, 1, mid);
	const int2 dim = make_int2(Kn, Pn);
	cuda_safe_call(__FILE__, __LINE__, [&] {
		cuda_mem_scope<float>(elastic_icdf_dev_p, pitch, dim, [&](float** elastic_p) {
		cuda_mem_scope<float>(inelastic_icdf_dev_p, pitch, dim, [&](float** inelastic_p) {
			float** K = new float*[dim.y];
			float** P = new float*[dim.y];
			K[0] = new float[dim.y*dim.x];
			P[0] = new float[dim.y*dim.x];
			for(int y = 0; y < dim.y; y++) {
				K[y] = K[0]+y*dim.x;
				P[y] = P[0]+y*dim.x;
				for(int x = 0; x < dim.x; x++) {
					K[y][x] = K1*std::exp(1.0*x/(dim.x-1)*std::log(K2/K1));
					P[y][x] = 1.0*y/(dim.y-1);
				}
			}
			callback(const_cast<const float**>(K), const_cast<const float**>(P), elastic_p, inelastic_p, dim);
			delete[] P[0];
			delete[] P;
			delete[] K[0];
			delete[] K;
		});
		});
	});
}

__host__ void material_struct::set_ionization(int mid, std::function<void(const float** K, const float** P, float** binding_p, int2 dim)> callback) {
	float* binding_dev_p = cuda_make_ptr<float>(ionization_dev_p, pitch, Pn+1, 1, mid);
	const int2 dim = make_int2(Kn, Pn);
	cuda_safe_call(__FILE__, __LINE__, [&] {
		cuda_mem_scope<float>(binding_dev_p, pitch, dim, [&](float** binding_p) {
			float** K = new float*[dim.y];
			float** P = new float*[dim.y];
			K[0] = new float[dim.y*dim.x];
			P[0] = new float[dim.y*dim.x];
			for(int y = 0; y < dim.y; y++) {
				K[y] = K[0]+y*dim.x;
				P[y] = P[0]+y*dim.x;
				for(int x = 0; x < dim.x; x++) {
					K[y][x] = K1*std::exp(1.0*x/(dim.x-1)*std::log(K2/K1));
					P[y][x] = 1.0*y/(dim.y-1);
				}
			}
			callback(const_cast<const float**>(K), const_cast<const float**>(P), binding_p, dim);
			delete[] P[0];
			delete[] P;
			delete[] K[0];
			delete[] K;
		});
	});
}