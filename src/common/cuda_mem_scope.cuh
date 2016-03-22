/**
 * @file src/common/cuda_mem_scope.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__CUDA_MEM_SCOPE__HEADER_INCLUDED
#define eSCATTER__COMMON__CUDA_MEM_SCOPE__HEADER_INCLUDED

#include <array>
#include <functional>

template<typename T>
__host__ void cuda_mem_scope(T* dev_p, int dim, std::function<void(T* host_p)>);
template<typename T>
__host__ void cuda_mem_scope(T* dev_p, int pitch, int2 dim, std::function<void(T** host_p)>);
template<typename T>
__host__ void cuda_mem_scope(T* dev_p, int pitch, int height, int3 dim, std::function<void(T*** host_p)>);

#include <common/cuda_mem_scope.cuinc>

#endif
