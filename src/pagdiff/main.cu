/**
 * @file src/pagdiff/pagdiff.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <array>
#include <cfloat>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct pagdiff_arg_t {
    float* grid;
    float3 grid_min;
    float3 grid_max;
    int3 grid_dim;
    float diff_radius;
    float thickness;
    float* pag_x;
    float* pag_y;
    float* pag_z;
    int pag_count;
};

__global__ void pagdiff_kernel(pagdiff_arg_t arg) {
    const int i = threadIdx.x+blockDim.x*threadIdx.y+blockDim.x*blockDim.y*threadIdx.z;
    const int ni = blockDim.x*blockDim.y*blockDim.z;
    const int jx = threadIdx.x+blockIdx.x*blockDim.x;
    const int jy = threadIdx.y+blockIdx.y*blockDim.y;
    const int jz = threadIdx.z+blockIdx.z*blockDim.z;
    bool in_grid = true;
    if((jx >= arg.grid_dim.x) || (jy >= arg.grid_dim.y) || (jz >= arg.grid_dim.z))
        in_grid = false;
    float sx, sy, sz;
    float x, y, z, w;
    if(in_grid) {
        sx = (1.0f*jx)/(arg.grid_dim.x-1);
        sy = (1.0f*jy)/(arg.grid_dim.y-1);
        sz = (1.0f*jz)/(arg.grid_dim.z-1);
        x = (1.0f-sx)*arg.grid_min.x+sx*arg.grid_max.x;
        y = (1.0f-sy)*arg.grid_min.y+sy*arg.grid_max.y;
        z = (1.0f-sz)*arg.grid_min.z+sz*arg.grid_max.z;
        w = -0.5f/(arg.diff_radius*arg.diff_radius);
    }
    float amplitude;
    if(in_grid)
        amplitude = arg.grid[jx+arg.grid_dim.x*jy+arg.grid_dim.x*arg.grid_dim.y*jz];
    __shared__ float _xi[1024];
    __shared__ float _yi[1024];
    __shared__ float _zi[1024];
    for(int j = 0; j < arg.pag_count; j += ni) {
        __syncthreads();
        if(i+j < arg.pag_count) {
            _xi[i] = arg.pag_x[i+j];
            _yi[i] = arg.pag_y[i+j];
            _zi[i] = arg.pag_z[i+j];
        }
        __syncthreads();
        if(in_grid) {
            const int nk = min(ni, arg.pag_count-j);
            for(int k = 0; k < nk; k++) {
                const float xi = _xi[k];
                const float yi = _yi[k];
                const float zi = _zi[k];
                const float dx = x-xi;
                const float dy = y-yi;
                const float dr2 = dx*dx+dy*dy;
                float dz;
                // direct amplitude
                dz = z-zi;
                amplitude += __expf(w*(dr2+dz*dz));
                // substrate reflection amplitude
                dz = z+zi;
                amplitude += __expf(w*(dr2+dz*dz));
                // vacuum reflection amplitude
                dz = z+zi-2.0f*arg.thickness;
                amplitude += __expf(w*(dr2+dz*dz));
            }
        }
    }
    if(in_grid)
        arg.grid[jx+arg.grid_dim.x*jy+arg.grid_dim.x*arg.grid_dim.y*jz] = amplitude;
}

__host__ void cuda_safe_call(std::function<void(void)> callback) {
    try {
        callback();
    } catch(...) {
        cuda_safe_call([&](){});
        throw;
    }
    if(cudaPeekAtLastError() != cudaSuccess) {
        std::ostringstream oss;
        oss << cudaGetErrorString(cudaGetLastError());
        throw std::runtime_error(oss.str());
    }
}

int main(const int argc, char* argv[]) {
    pagdiff_arg_t arg;
    arg.grid_min = make_float3(-64, -400, 0);
    arg.grid_max = make_float3(64, 400, 100);
    arg.grid_dim = make_int3(256, 1024, 128);
    arg.diff_radius = 5;
    arg.thickness = 100;
    const dim3 block_size = dim3(8, 8, 8);
    const dim3 block_count = dim3(
        1+arg.grid_dim.x/block_size.x,
        1+arg.grid_dim.y/block_size.y,
        1+arg.grid_dim.z/block_size.z
    );
    const int batch_size = 10*block_size.x*block_size.y*block_size.z;
    cuda_safe_call([&]{
        const size_t n = arg.grid_dim.x*arg.grid_dim.y*arg.grid_dim.z;
        cudaMalloc(&arg.grid, n*sizeof(float));
        cudaMemset(arg.grid, 0, n*sizeof(float));
        cudaMalloc(&arg.pag_x, batch_size*sizeof(float));
        cudaMalloc(&arg.pag_y, batch_size*sizeof(float));
        cudaMalloc(&arg.pag_z, batch_size*sizeof(float));
    });
    std::vector<float> pag_x(batch_size);
    std::vector<float> pag_y(batch_size);
    std::vector<float> pag_z(batch_size);
    size_t count = 0;
    while(!std::cin.eof()) {
        std::clog << "\r >> computing " << count << std::flush;
        arg.pag_count = 0;
        while(arg.pag_count < batch_size) {
            std::string dummy;
            std::cin >> dummy >> dummy >> dummy;
            if(std::cin.eof())
                break;
            std::string rx, ry, rz;
            std::cin >> rx >> ry >> rz;
            pag_x[arg.pag_count] = std::stof(rx);
            pag_y[arg.pag_count] = std::stof(ry);
            pag_z[arg.pag_count] = std::stof(rz);
            if((pag_x[arg.pag_count] >= arg.grid_min.x) && (pag_x[arg.pag_count] <= arg.grid_max.x))
            if((pag_y[arg.pag_count] >= arg.grid_min.y) && (pag_y[arg.pag_count] <= arg.grid_max.y))
            if((pag_z[arg.pag_count] >= arg.grid_min.z) && (pag_z[arg.pag_count] <= arg.grid_max.z))
                arg.pag_count++;
        }
        cuda_safe_call([&]{
            cudaMemcpy(arg.pag_x, pag_x.data(), arg.pag_count*sizeof(float), cudaMemcpyHostToDevice);
            cudaMemcpy(arg.pag_y, pag_y.data(), arg.pag_count*sizeof(float), cudaMemcpyHostToDevice);
            cudaMemcpy(arg.pag_z, pag_z.data(), arg.pag_count*sizeof(float), cudaMemcpyHostToDevice);
            pagdiff_kernel<<<block_count,block_size>>>(arg);
        });
        count += arg.pag_count;
    }
    std::clog << std::endl;
    std::vector<float> grid(arg.grid_dim.x*arg.grid_dim.y*arg.grid_dim.z);
    cuda_safe_call([&]{
        cudaMemcpy(grid.data(), arg.grid, grid.size()*sizeof(float), cudaMemcpyDeviceToHost);
    });
    std::cout.precision(FLT_DIG);
    std::cout << std::scientific;
    for(int iz = 0; iz < arg.grid_dim.z; iz++)
    for(int iy = 0; iy < arg.grid_dim.y; iy++) {
        for(int ix = 0; ix < arg.grid_dim.x; ix++) {
            if(ix > 0)
                std::cout << " ";
            std::cout << grid[ix+arg.grid_dim.x*iy+arg.grid_dim.x*arg.grid_dim.y*iz];
        }
        std::cout << std::endl;
    }
    cuda_safe_call([&]{
        cudaFree(arg.grid);
        cudaFree(arg.pag_x);
        cudaFree(arg.pag_y);
        cudaFree(arg.pag_z);
    });
    return EXIT_SUCCESS;
}