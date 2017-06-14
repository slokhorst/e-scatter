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
#include <csread/units/unit_system.h>

__host__ cuda_material_struct cuda_material_struct::create(int capacity) {
    cuda_material_struct mstruct;
    mstruct.capacity = capacity;
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        size_t _pitch;
        cudaMalloc(&mstruct.fermi_dev_p, capacity*sizeof(float));
        cudaMalloc(&mstruct.barrier_dev_p, capacity*sizeof(float));
        cudaMalloc(&mstruct.band_gap_dev_p, capacity*sizeof(float));
        cudaMalloc(&mstruct.band_edge_dev_p, capacity*sizeof(float));
        cudaMalloc(&mstruct.effective_mass_dev_p, capacity*sizeof(float));
        cudaMalloc(&mstruct.phonon_loss_dev_p, capacity*sizeof(float));
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
        cudaFree(mstruct.band_gap_dev_p);
        cudaFree(mstruct.band_edge_dev_p);
        cudaFree(mstruct.effective_mass_dev_p);
        cudaFree(mstruct.phonon_loss_dev_p);
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
            fermi_p[i] = (_material.get_fermi()/units::eV).value;
        });
        cuda_mem_scope<float>(barrier_dev_p, capacity, [&](float* barrier_p) {
            barrier_p[i] = (_material.get_barrier()/units::eV).value;
        });
        cuda_mem_scope<float>(band_gap_dev_p, capacity, [&](float* band_gap_p) {
            if(_material.get_conductor_type() != material::CND_METAL)
                band_gap_p[i] = (_material.get_band_gap()/units::eV).value;
            else
                band_gap_p[i] = -1;
        });
        cuda_mem_scope<float>(band_edge_dev_p, capacity, [&](float* band_edge_p) {
            #warning "band edge is determined for silicon and pmma by using string comparison"
            if(_material.get_name() == "silicon" || _material.get_name() == "silicon-surface")
                band_edge_p[i] = (_material.get_barrier()/units::eV).value - 4.05f;
            else if(_material.get_name() == "pmma" || _material.get_name() == "pmma-surface")
                band_edge_p[i] = (_material.get_barrier() / units::eV).value - 2.5f;
            else
                band_edge_p[i] = 0;
        });
        cuda_mem_scope<float>(effective_mass_dev_p, capacity, [&](float* effective_mass_p) {
            #warning "effective mass is determined for silicon and pmma by using string comparison"
            if(_material.get_name() == "silicon" || _material.get_name() == "silicon-surface")
                effective_mass_p[i] = 1.09f;
            else if(_material.get_name() == "pmma" || _material.get_name() == "pmma-surface")
                effective_mass_p[i] = 1.0f;
            else
                effective_mass_p[i] = 1.0f;
        });
        cuda_mem_scope<float>(phonon_loss_dev_p, capacity, [&](float* phonon_loss_p) {
            phonon_loss_p[i] = (_material.get_phonon_loss()/units::eV).value;
        });
    });

    cuda_safe_call(__FILE__, __LINE__, [&]() {
        const auto elastic_imfp_tbl = _material.get_elastic_imfp(K_min, K_max, K_cnt);
        const auto inelastic_imfp_tbl = _material.get_inelastic_imfp(K_min, K_max, K_cnt);

        float* elastic_imfp_dev_p = cuda_make_ptr<float>(elastic_dev_p, pitch, P_cnt+1, 0, i);
        float* inelastic_imfp_dev_p = cuda_make_ptr<float>(inelastic_dev_p, pitch, P_cnt+1, 0, i);
        cuda_mem_scope<float>(elastic_imfp_dev_p, K_cnt, [&](float* elastic_imfp_p) {
        cuda_mem_scope<float>(inelastic_imfp_dev_p, K_cnt, [&](float* inelastic_imfp_p) {
            for(int x = 0; x < K_cnt; x++) {
                elastic_imfp_p[x] = elastic_imfp_tbl(x);
                inelastic_imfp_p[x] = inelastic_imfp_tbl(x);
            }
        });
        });
    });
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        const auto elastic_icdf_tbl = _material.get_elastic_angle_icdf(K_min, K_max, K_cnt, P_cnt);
        const auto inelastic_icdf_tbl = _material.get_inelastic_w0_icdf(K_min, K_max, K_cnt, P_cnt);

        float* elastic_icdf_dev_p = cuda_make_ptr<float>(elastic_dev_p, pitch, P_cnt+1, 1, i);
        float* inelastic_icdf_dev_p = cuda_make_ptr<float>(inelastic_dev_p, pitch, P_cnt+1, 1, i);
        cuda_mem_scope<float>(elastic_icdf_dev_p, pitch, make_int2(K_cnt, P_cnt), [&](float** elastic_icdf_p) {
        cuda_mem_scope<float>(inelastic_icdf_dev_p, pitch, make_int2(K_cnt, P_cnt), [&](float** inelastic_icdf_p) {
            for(int y = 0; y < P_cnt; y++) {
                for(int x = 0; x < K_cnt; x++) {
                    const double K = inelastic_icdf_tbl.get_x(x);
                    elastic_icdf_p[y][x] = std::cos(max(0.0, min(M_PI, elastic_icdf_tbl(x, y))));
                    inelastic_icdf_p[y][x] = std::log(max(0.0, min(K-(_material.get_fermi()/units::eV).value, inelastic_icdf_tbl(x, y))));
                }
            }
        });
        });
    });
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        const auto binding_tbl = _material.get_ionization_icdf(K_min, K_max, K_cnt, P_cnt);

        float* binding_dev_p = cuda_make_ptr<float>(ionization_dev_p, pitch, P_cnt+1, 0, i);
        cuda_mem_scope<float>(binding_dev_p, pitch, make_int2(K_cnt, P_cnt), [&](float** binding_p) {
            for(int y = 0; y < P_cnt; y++) {
                for(int x = 0; x < K_cnt; x++) {
                    binding_p[y][x] = binding_tbl(x, y);
                }
            }
        });
    });
}