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

    int* octree_dev_p;
    int octree_pitch;
    float3 AABB_center;
    float3 AABB_halfsize;
    int occupancy;

    int* material_idx_in_dev_p; int* material_idx_out_dev_p;
    float* triangle_r0x_dev_p; float* triangle_r0y_dev_p; float* triangle_r0z_dev_p;
    float* triangle_e1x_dev_p; float* triangle_e1y_dev_p; float* triangle_e1z_dev_p;
    float* triangle_e2x_dev_p; float* triangle_e2y_dev_p; float* triangle_e2z_dev_p;

private:
    __host__ __device__ cuda_geometry_struct() = default;
};

#endif