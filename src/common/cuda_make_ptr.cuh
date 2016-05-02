/**
 * @file src/common/cuda_make_ptr.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__CUDA_MAKE_PTR__HEADER_INCLUDED
#define eSCATTER__COMMON__CUDA_MAKE_PTR__HEADER_INCLUDED

template<typename T>
__host__ __device__ const T* cuda_make_ptr(const T* ptr, int pitch, int iy);
template<typename T>
__host__ __device__ const T* cuda_make_ptr(const T* ptr, int pitch, int height, int iy, int iz);

template<typename T>
__host__ __device__ T* cuda_make_ptr(T* ptr, int pitch, int iy, int iz);
template<typename T>
__host__ __device__ T* cuda_make_ptr(T* ptr, int pitch, int height, int iy, int iz);

#include <common/cuda_make_ptr.cuinc>

#endif
