/**
 * @file src/cdsem/cuda_material_struct.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "cuda_material_struct.cuh"
#include <cuda_common/cuda_make_ptr.cuh>
#include <cuda_common/cuda_mem_scope.cuh>
#include <cuda_common/cuda_safe_call.cuh>
#include <common/constant.hh>

__host__ cuda_material_struct cuda_material_struct::create(int capacity) {
    cuda_material_struct mstruct;
    mstruct.capacity = capacity;
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        size_t _pitch;
        cudaMalloc(&mstruct.fermi_dev_p, capacity*sizeof(float));
        cudaMalloc(&mstruct.barrier_dev_p, capacity*sizeof(float));
        cudaMalloc(&mstruct.bandgap_dev_p, capacity*sizeof(float));
        cudaMalloc(&mstruct.phononloss_dev_p, capacity*sizeof(float));
        cudaMallocPitch(&mstruct.elastic_dev_p, &_pitch, mstruct.K_cnt*sizeof(float), (mstruct.P_cnt+1)*capacity);
        cudaMallocPitch(&mstruct.inelastic_dev_p, &_pitch, mstruct.K_cnt*sizeof(float), (mstruct.P_cnt+1)*capacity);
        cudaMallocPitch(&mstruct.ionization_dev_p, &_pitch, mstruct.K_cnt*sizeof(float), (mstruct.P_cnt+1)*capacity);
        mstruct.pitch = _pitch;
    });
    return mstruct;
};

__host__ void cuda_material_struct::release(cuda_material_struct& mstruct) {
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaFree(mstruct.fermi_dev_p);
        cudaFree(mstruct.barrier_dev_p);
        cudaFree(mstruct.bandgap_dev_p);
        cudaFree(mstruct.phononloss_dev_p);
        cudaFree(mstruct.elastic_dev_p);
        cudaFree(mstruct.inelastic_dev_p);
        cudaFree(mstruct.ionization_dev_p);
    });
}

__host__ void cuda_material_struct::assign(int i, const material& _material) {
    if((i < 0) || (i >= capacity))
        return;
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cuda_mem_scope<float>(fermi_dev_p, capacity, [&](float* fermi_p) {
            fermi_p[i] = _material.fermi()/constant::ec;
        });
        cuda_mem_scope<float>(barrier_dev_p, capacity, [&](float* barrier_p) {
            barrier_p[i] = _material.barrier()/constant::ec;
        });
        cuda_mem_scope<float>(bandgap_dev_p, capacity, [&](float* bandgap_p) {
            if(_material.bandgap().is_defined())
                bandgap_p[i] = _material.bandgap()()/constant::ec;
            else
                bandgap_p[i] = -1;
        });
        cuda_mem_scope<float>(phononloss_dev_p, capacity, [&](float* phononloss_p) {
            phononloss_p[i] = _material.phononloss()/constant::ec;
        });
    });
    auto __logspace_K_at = [&](int x) {
        return K_min*std::exp(1.0*x/(K_cnt-1)*std::log(K_max/K_min));
    };
    auto __linspace_P_at = [&](int y) {
        return 1.0*y/(P_cnt-1);
    };
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        float* elastic_imfp_dev_p = cuda_make_ptr<float>(elastic_dev_p, pitch, P_cnt+1, 0, i);
        float* inelastic_imfp_dev_p = cuda_make_ptr<float>(inelastic_dev_p, pitch, P_cnt+1, 0, i);
        cuda_mem_scope<float>(elastic_imfp_dev_p, K_cnt, [&](float* elastic_imfp_p) {
        cuda_mem_scope<float>(inelastic_imfp_dev_p, K_cnt, [&](float* inelastic_imfp_p) {
            for(int x = 0; x < K_cnt; x++) {
                const double K = __logspace_K_at(x)*constant::ec;
                elastic_imfp_p[x] = std::log(_material.density()*_material.elastic_tcs(K)*1e-9);
                inelastic_imfp_p[x] = std::log(_material.density()*_material.inelastic_tcs(K)*1e-9);
            }
        });
        });
    });
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        float* elastic_icdf_dev_p = cuda_make_ptr<float>(elastic_dev_p, pitch, P_cnt+1, 1, i);
        float* inelastic_icdf_dev_p = cuda_make_ptr<float>(inelastic_dev_p, pitch, P_cnt+1, 1, i);
        cuda_mem_scope<float>(elastic_icdf_dev_p, pitch, make_int2(K_cnt, P_cnt), [&](float** elastic_icdf_p) {
        cuda_mem_scope<float>(inelastic_icdf_dev_p, pitch, make_int2(K_cnt, P_cnt), [&](float** inelastic_icdf_p) {
            for(int y = 0; y < P_cnt; y++) {
                const double P =  __linspace_P_at(y);
                for(int x = 0; x < K_cnt; x++) {
                    const double K = __logspace_K_at(x)*constant::ec;
                    elastic_icdf_p[y][x] = std::cos(_material.elastic_icdf(K, P));
                    inelastic_icdf_p[y][x] = std::log(_material.inelastic_icdf(K, P)/constant::ec);
                }
            }
        });
        });
    });
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        float* binding_dev_p = cuda_make_ptr<float>(ionization_dev_p, pitch, P_cnt+1, 0, i);
        cuda_mem_scope<float>(binding_dev_p, pitch, make_int2(K_cnt, P_cnt), [&](float** binding_p) {
            for(int y = 0; y < P_cnt; y++) {
                const double P = __linspace_P_at(y);
                for(int x = 0; x < K_cnt; x++) {
                    const double omega0 = __logspace_K_at(x);
                    const double margin = 10; // magic number in accordance with Kieft & Bosch code
                    double binding = -1;
                    if(omega0 > 100) {
                        binding = _material.ionization_energy((omega0+margin)*constant::ec, P)/constant::ec;
                        if(binding < 50)
                            binding = -1;
                    }
                    if(binding < 0) {
                        binding = _material.outer_shell_ionization_energy(omega0*constant::ec)/constant::ec;
                        if(binding < 0)
                            binding = -1;
                    }
                    binding_p[y][x] = binding;
                }
            }
        });
    });
}