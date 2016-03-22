/**
 * @file src/common/cuda_vec3_math.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__CUDA_VEC3_MATH__HEADER_INCLUDED
#define eSCATTER__COMMON__CUDA_VEC3_MATH__HEADER_INCLUDED

inline __host__ __device__ float3 operator+(float3 a, float3 b);
inline __host__ __device__ float3 operator-(float3 a, float3 b);
inline __host__ __device__ float3 operator*(float3 a, float b);
inline __host__ __device__ float3 operator*(float a, float3 b);
inline __host__ __device__ float3 operator-(float3 a);
inline __host__ __device__ float3 cross_product(float3, float3);
inline __host__ __device__ float dot_product(float3, float3);
inline __host__ __device__ float norm2(float3);
inline __host__ __device__ float3 normalize(float3);

#include <common/cuda_vec3_math.cuinc>

#endif