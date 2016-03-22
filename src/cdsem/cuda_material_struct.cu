/**
 * @file src/cdsem/cuda_material_struct.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "cuda_material_struct.cuh"
#include <common/cuda_make_ptr.cuh>
#include <common/cuda_mem_scope.cuh>
#include <common/cuda_safe_call.cuh>
#include <common/constant.h>

__host__ cuda_material_struct cuda_material_struct::create(int n) {
	cuda_material_struct mstruct;
	mstruct.n = n;
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		size_t _pitch;
		cudaMalloc(&mstruct.fermi_dev_p, n*sizeof(float));
		cudaMalloc(&mstruct.barrier_dev_p, n*sizeof(float));
		cudaMalloc(&mstruct.bandgap_dev_p, n*sizeof(float));
		cudaMallocPitch(&mstruct.elastic_dev_p, &_pitch, mstruct.Kn*sizeof(float), (mstruct.Pn+1)*n);
		cudaMallocPitch(&mstruct.inelastic_dev_p, &_pitch, mstruct.Kn*sizeof(float), (mstruct.Pn+1)*n);
		cudaMallocPitch(&mstruct.ionization_dev_p, &_pitch, mstruct.Kn*sizeof(float), mstruct.Pn*n);
		mstruct.pitch = _pitch;
	});
	return mstruct;
};

__host__ void cuda_material_struct::release(cuda_material_struct& mstruct) {
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cudaFree(mstruct.fermi_dev_p);
		cudaFree(mstruct.barrier_dev_p);
		cudaFree(mstruct.bandgap_dev_p);
		cudaFree(mstruct.elastic_dev_p);
		cudaFree(mstruct.inelastic_dev_p);
		cudaFree(mstruct.ionization_dev_p);
	});
}

__host__ void cuda_material_struct::assign(int i, material& mat) {
	if((i < 0) || (i >= n))
		return;
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		cuda_mem_scope<float>(fermi_dev_p, n, [&](float* fermi_p) {
			fermi_p[i] = mat.fermi()/constant::ec;
		});
		cuda_mem_scope<float>(barrier_dev_p, n, [&](float* barrier_p) {
			barrier_p[i] = mat.barrier()/constant::ec;
		});
		cuda_mem_scope<float>(bandgap_dev_p, n, [&](float* bandgap_p) {
			if(mat.bandgap().is_defined())
				bandgap_p[i] = mat.bandgap()()/constant::ec;
			else
				bandgap_p[i] = -1;
		});
	});
	auto __K_at = [&](int x) {
		return K1*std::exp(1.0*x/(Kn-1)*std::log(K2/K1));
	};
	auto __P_at = [&](int y) {
		return 1.0*y/(Pn-1);
	};
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		float* elastic_imfp_dev_p = cuda_make_ptr<float>(elastic_dev_p, pitch, Pn+1, 0, i);
		float* inelastic_imfp_dev_p = cuda_make_ptr<float>(inelastic_dev_p, pitch, Pn+1, 0, i);
		cuda_mem_scope<float>(elastic_imfp_dev_p, Kn, [&](float* elastic_imfp_p) {
		cuda_mem_scope<float>(inelastic_imfp_dev_p, Kn, [&](float* inelastic_imfp_p) {
			for(int x = 0; x < Kn; x++) {
				elastic_imfp_p[x] = std::log(mat.density()*mat.elastic_tcs(__K_at(x)*constant::ec)*1e-9);
				inelastic_imfp_p[x] = std::log(mat.density()*mat.inelastic_tcs(__K_at(x)*constant::ec)*1e-9);
			}
		});
		});
	});
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		float* elastic_icdf_dev_p = cuda_make_ptr<float>(elastic_dev_p, pitch, Pn+1, 1, i);
		float* inelastic_icdf_dev_p = cuda_make_ptr<float>(inelastic_dev_p, pitch, Pn+1, 1, i);
		cuda_mem_scope<float>(elastic_icdf_dev_p, pitch, make_int2(Kn, Pn), [&](float** elastic_icdf_p) {
		cuda_mem_scope<float>(inelastic_icdf_dev_p, pitch, make_int2(Kn, Pn), [&](float** inelastic_icdf_p) {
			for(int y = 0; y < Pn; y++)
			for(int x = 0; x < Kn; x++) {
					elastic_icdf_p[y][x] = std::cos(mat.elastic_dcs(__K_at(x)*constant::ec, __P_at(y)));
					inelastic_icdf_p[y][x] = mat.inelastic_dcs(__K_at(x)*constant::ec, __P_at(y))/constant::ec;
			}
		});
		});
	});
	cuda_safe_call(__FILE__, __LINE__, [&]() {
		float* binding_dev_p = cuda_make_ptr<float>(ionization_dev_p, pitch, Kn, 0, i);
		cuda_mem_scope<float>(binding_dev_p, pitch, make_int2(Kn, Pn), [&](float** binding_p) {
			for(int y = 0; y < Pn; y++)
			for(int x = 0; x < Kn; x++) {
				binding_p[y][x] = mat.ionization_energy(__K_at(x)*constant::ec, __P_at(y))/constant::ec;
			}
		});
	});
}