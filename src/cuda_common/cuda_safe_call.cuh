/**
 * @file src/common/cuda_safe_call.cuh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__CUDA_SAFE_CALL__HEADER_INCLUDED
#define eSCATTER__COMMON__CUDA_SAFE_CALL__HEADER_INCLUDED

#include <functional>

__host__ void cuda_safe_call(const char* file, int line, std::function<void(void)>);

#endif
