/**
 * @file src/common/cuda_safe_call.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <common/cuda_safe_call.cuh>
#include <sstream>
#include <stdexcept>

__host__ void cuda_safe_call(const char* file, int line, std::function<void(void)> callback) {
    try {
        callback();
    } catch(...) {
        cuda_safe_call(file, line, [&](){});
        throw;
    }
    if(cudaPeekAtLastError() != cudaSuccess) {
        std::ostringstream oss;
        oss << cudaGetErrorString(cudaGetLastError());
        oss << " in '" << file << "' at line " << line;
        throw std::runtime_error(oss.str());
    }
}
