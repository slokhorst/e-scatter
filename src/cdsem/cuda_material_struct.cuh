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
    __host__ static cuda_material_struct create(int n);
    __host__ static void release(cuda_material_struct&);
    __host__ void assign(int i, material& mat);
    int n;
    const int Kn = 64;
    const int Pn = 64;
    const float K1 = 1;
    const float K2 = 10e3;
    float* fermi_dev_p;
    float* barrier_dev_p;
    float* bandgap_dev_p;
    float* elastic_dev_p;
    float* inelastic_dev_p;
    float* ionization_dev_p;
    int pitch;
private:
    __host__ __device__ cuda_material_struct() = default;
};

#endif