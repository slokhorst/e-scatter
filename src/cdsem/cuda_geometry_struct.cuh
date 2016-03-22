/**
 * @file src/cdsem/cuda_geometry_struct.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__CUDA_GEOMETRY_STRUCT__HEADER_INCLUDED
#define eSCATTER__CDSEM__CUDA_GEOMETRY_STRUCT__HEADER_INCLUDED

#include "trimesh.h"

struct cuda_geometry_struct {
	__host__ static cuda_geometry_struct create(const trimesh&, int cell_count);
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
	float3 org;
	int4 dim;
	float3 cell;
	int* map_dev_p;
	int* in_dev_p; int* out_dev_p;
	float* Ax_dev_p; float* Ay_dev_p; float* Az_dev_p;
	float* Bx_dev_p; float* By_dev_p; float* Bz_dev_p;
	float* Cx_dev_p; float* Cy_dev_p; float* Cz_dev_p;
	int pitch;
private:
	__host__ __device__ cuda_geometry_struct() = default;
};

#endif