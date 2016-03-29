/**
 * @file src/cdsem/cuda_geometry_struct.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "cuda_geometry_struct.cuh"
#include <algorithm>
#include <cfloat>
#include <functional>
#include <map>
#include <vector>
#include <common/cuda_mem_scope.cuh>
#include <common/cuda_safe_call.cuh>
#include "trigrid.hh"

__host__ cuda_geometry_struct cuda_geometry_struct::create(const trimesh& triangle_mesh, int cell_count) {
    cuda_geometry_struct gstruct;
    gstruct.org = make_float3(0, 0, 0);
    gstruct.dim = make_int4(0, 0, 0, 0);
    gstruct.cell = make_float3(0, 0, 0);
    gstruct.map_dev_p = nullptr;
    gstruct.in_dev_p = gstruct.out_dev_p = nullptr;
    gstruct.Ax_dev_p = gstruct.Ay_dev_p = gstruct.Az_dev_p = nullptr;
    gstruct.Bx_dev_p = gstruct.By_dev_p = gstruct.Bz_dev_p = nullptr;
    gstruct.Cx_dev_p = gstruct.Cy_dev_p = gstruct.Cz_dev_p = nullptr;
    gstruct.pitch = 0;
    if(triangle_mesh.empty())
        return gstruct;
    point3 grid_min = {0, 0, 0};
    point3 grid_max = {0, 0, 0};
    auto cit = triangle_mesh.cbegin();
    std::tie(grid_min.x, grid_max.x) = std::minmax({cit->A.x, cit->B.x, cit->C.x});
    std::tie(grid_min.y, grid_max.y) = std::minmax({cit->A.y, cit->B.y, cit->C.y});
    std::tie(grid_min.z, grid_max.z) = std::minmax({cit->A.z, cit->B.z, cit->C.z});
    while(++cit != triangle_mesh.cend()) {
        std::tie(grid_min.x, grid_max.x) = std::minmax({grid_min.x, grid_max.x, cit->A.x, cit->B.x, cit->C.x});
        std::tie(grid_min.y, grid_max.y) = std::minmax({grid_min.y, grid_max.y, cit->A.y, cit->B.y, cit->C.y});
        std::tie(grid_min.z, grid_max.z) = std::minmax({grid_min.z, grid_max.z, cit->A.z, cit->B.z, cit->C.z});
    }
    grid_min.x -= 1; grid_min.y -= 1; grid_min.z -= 1;
    grid_max.x += 1; grid_max.y += 1; grid_max.z += 1;
    const point3 grid_size(grid_max.x-grid_min.x, grid_max.y-grid_min.y, grid_max.z-grid_min.z);
    const double grid_delta = std::pow(grid_size.x*grid_size.y*grid_size.z/cell_count, 1.0/3);
    const index3 grid_dim(std::ceil(grid_size.x/grid_delta), std::ceil(grid_size.y/grid_delta), std::ceil(grid_size.z/grid_delta));
    trigrid triangle_grid(grid_min, grid_max, grid_dim);
    for(auto cit = triangle_mesh.cbegin(); cit != triangle_mesh.cend(); cit++)
        triangle_grid.insert(*cit);
    gstruct.org = make_float3(grid_min.x, grid_min.y, grid_min.z);
    gstruct.dim = make_int4(grid_dim.x, grid_dim.y, grid_dim.z, 0);
    gstruct.cell = make_float3(grid_size.x/grid_dim.x, grid_size.y/grid_dim.y, grid_size.z/grid_dim.z);
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&gstruct.map_dev_p, gstruct.dim.x*gstruct.dim.y*gstruct.dim.z*sizeof(int));
        cuda_mem_scope<int>(gstruct.map_dev_p, gstruct.dim.x*gstruct.dim.y*gstruct.dim.z, [&](int* map_p) {
            for(int i = 0; i < gstruct.dim.x*gstruct.dim.y*gstruct.dim.z; i++)
                map_p[i] = -1;
        });
    });
    int2 size = make_int2(0, 0);
    triangle_grid.for_each_cell([&](const index3& i, const trigrid::triangle_p_vector& triangle_p_vec) {
        size.x = std::max(size.x, static_cast<int>(triangle_p_vec.size()));
        size.y++;
    });
    if(size.y == 0)
        return gstruct;
    gstruct.dim.w = size.x;
    const int _sizeof = std::max(sizeof(int), sizeof(float));
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        size_t _pitch;
        for(int** p : {&(gstruct.in_dev_p), &(gstruct.out_dev_p)})
            cudaMallocPitch(p, &_pitch, size.x*_sizeof, size.y);
        for(float** p : {&(gstruct.Ax_dev_p), &(gstruct.Ay_dev_p), &(gstruct.Az_dev_p)})
            cudaMallocPitch(p, &_pitch, size.x*_sizeof, size.y);
        for(float** p : {&(gstruct.Bx_dev_p), &(gstruct.By_dev_p), &(gstruct.Bz_dev_p)})
            cudaMallocPitch(p, &_pitch, size.x*_sizeof, size.y);
        for(float** p : {&(gstruct.Cx_dev_p), &(gstruct.Cy_dev_p), &(gstruct.Cz_dev_p)})
            cudaMallocPitch(p, &_pitch, size.x*_sizeof, size.y);
        gstruct.pitch = _pitch;
    });
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cuda_mem_scope<int>(gstruct.in_dev_p, gstruct.pitch*size.y/_sizeof, [&](int* in_p) {
        cuda_mem_scope<int>(gstruct.out_dev_p, gstruct.pitch*size.y/_sizeof, [&](int* out_p) {
            for(int i = 0; i < gstruct.pitch*size.y/_sizeof; i++)
                in_p[i] = out_p[i] = NOP;
        });
        });
        cuda_mem_scope<int>(gstruct.map_dev_p, gstruct.dim.x*gstruct.dim.y*gstruct.dim.z, [&](int* map_p) {
            cuda_mem_scope<int>(gstruct.in_dev_p, gstruct.pitch, size, [&](int** in_p) {
            cuda_mem_scope<int>(gstruct.out_dev_p, gstruct.pitch, size, [&](int** out_p) {
                cuda_mem_scope<float>(gstruct.Ax_dev_p, gstruct.pitch, size, [&](float** Ax_p) {
                cuda_mem_scope<float>(gstruct.Ay_dev_p, gstruct.pitch, size, [&](float** Ay_p) {
                cuda_mem_scope<float>(gstruct.Az_dev_p, gstruct.pitch, size, [&](float** Az_p) {
                    cuda_mem_scope<float>(gstruct.Bx_dev_p, gstruct.pitch, size, [&](float** Bx_p) {
                    cuda_mem_scope<float>(gstruct.By_dev_p, gstruct.pitch, size, [&](float** By_p) {
                    cuda_mem_scope<float>(gstruct.Bz_dev_p, gstruct.pitch, size, [&](float** Bz_p) {
                        cuda_mem_scope<float>(gstruct.Cx_dev_p, gstruct.pitch, size, [&](float** Cx_p) {
                        cuda_mem_scope<float>(gstruct.Cy_dev_p, gstruct.pitch, size, [&](float** Cy_p) {
                        cuda_mem_scope<float>(gstruct.Cz_dev_p, gstruct.pitch, size, [&](float** Cz_p) {
                            int y = 0;
                            triangle_grid.for_each_cell([&](const index3& i, const trigrid::triangle_p_vector& triangle_p_vec) {
                                for(size_t x = 0; x < triangle_p_vec.size(); x++) {
                                    in_p[y][x] = triangle_p_vec[x]->in;
                                    out_p[y][x] = triangle_p_vec[x]->out;
                                    Ax_p[y][x] = triangle_p_vec[x]->A.x;
                                    Ay_p[y][x] = triangle_p_vec[x]->A.y;
                                    Az_p[y][x] = triangle_p_vec[x]->A.z;
                                    Bx_p[y][x] = triangle_p_vec[x]->B.x;
                                    By_p[y][x] = triangle_p_vec[x]->B.y;
                                    Bz_p[y][x] = triangle_p_vec[x]->B.z;
                                    Cx_p[y][x] = triangle_p_vec[x]->C.x;
                                    Cy_p[y][x] = triangle_p_vec[x]->C.y;
                                    Cz_p[y][x] = triangle_p_vec[x]->C.z;
                                }
                                const int gid = i.x+gstruct.dim.x*i.y+gstruct.dim.x*gstruct.dim.y*i.z;
                                map_p[gid] = y++;
                            });
                        });
                        });
                        });
                    });
                    });
                    });
                });
                });
                });
            });
            });
        });
    });
    return gstruct;
}

__host__ void cuda_geometry_struct::release(cuda_geometry_struct& gstruct) {
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaFree(gstruct.map_dev_p);
    });
    gstruct.map_dev_p = nullptr;
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        for(int* p : {gstruct.in_dev_p, gstruct.out_dev_p})
            cudaFree(p);
        for(float* p : {gstruct.Ax_dev_p, gstruct.Ay_dev_p, gstruct.Az_dev_p})
            cudaFree(p);
        for(float* p : {gstruct.Bx_dev_p, gstruct.By_dev_p, gstruct.Bz_dev_p})
            cudaFree(p);
        for(float* p : {gstruct.Cx_dev_p, gstruct.Cy_dev_p, gstruct.Cz_dev_p})
            cudaFree(p);
    });
    gstruct.in_dev_p = gstruct.out_dev_p = nullptr;
    gstruct.Ax_dev_p = gstruct.Ay_dev_p = gstruct.Az_dev_p = nullptr;
    gstruct.Bx_dev_p = gstruct.By_dev_p = gstruct.Bz_dev_p = nullptr;
    gstruct.Cx_dev_p = gstruct.Cy_dev_p = gstruct.Cz_dev_p = nullptr;
}