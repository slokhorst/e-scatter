/**
 * @file src/cdsem/cuda_material_struct.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__CUDA_MATERIAL_STRUCT__HEADER_INCLUDED
#define eSCATTER__CDSEM__CUDA_MATERIAL_STRUCT__HEADER_INCLUDED

#include <functional>
#include "material.hh"

class cuda_material_struct {
public:
    __host__ static cuda_material_struct create(int capacity);
    __host__ static void release(cuda_material_struct&);

    __host__ void assign(int i, material& mat);

    int capacity;
    const float2 K_energy_range = make_float2(1, 10e3);
    const int2 table_dim = make_int2(64, 64); // (kinetic energy, cumulative probability)
    int table_pitch;
    float* fermi_dev_p;
    float* barrier_dev_p;
    float* bandgap_dev_p;
    float* elastic_dev_p;
    float* inelastic_dev_p;
    float* ionization_dev_p;

private:
    __host__ __device__ cuda_material_struct() = default;
};

#endif