/**
 * @file src/cdsem/cuda_geometry_struct.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__CUDA_GEOMETRY_STRUCT__HEADER_INCLUDED
#define eSCATTER__CDSEM__CUDA_GEOMETRY_STRUCT__HEADER_INCLUDED

#include "octree.hh"

struct cuda_geometry_struct {
    __host__ static cuda_geometry_struct create(const octree&);
    __host__ static void release(cuda_geometry_struct&);

    enum id_enum : int {
        NOP = -128,
        TERMINATOR = -127,
        DETECTOR = -126,
        DETECTOR_LT50 = -125,
        DETECTOR_GE50 = -124,
        VACUUM = -123,
        MIRROR = -122
    };

    int* octree_dev_p;
    int octree_pitch;
    float3 root_center;
    float3 root_size;
    int occupancy;

    int* material_idx_in_dev_p; int* material_idx_out_dev_p;
    float* triangle_Ax_dev_p; float* triangle_Ay_dev_p; float* triangle_Az_dev_p;
    float* triangle_Bx_dev_p; float* triangle_By_dev_p; float* triangle_Bz_dev_p;
    float* triangle_Cx_dev_p; float* triangle_Cy_dev_p; float* triangle_Cz_dev_p;

private:
    __host__ __device__ cuda_geometry_struct() = default;
};

#endif