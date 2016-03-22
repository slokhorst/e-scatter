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
    point3 min = {0, 0, 0};
    point3 max = {0, 0, 0};
    auto cit = triangle_mesh.cbegin();
    std::tie(min.x, max.x) = std::minmax({cit->A.x, cit->B.x, cit->C.x});
    std::tie(min.y, max.y) = std::minmax({cit->A.y, cit->B.y, cit->C.y});
    std::tie(min.z, max.z) = std::minmax({cit->A.z, cit->B.z, cit->C.z});
    while(++cit != triangle_mesh.cend()) {
        std::tie(min.x, max.x) = std::minmax({min.x, max.x, cit->A.x, cit->B.x, cit->C.x});
        std::tie(min.y, max.y) = std::minmax({min.y, max.y, cit->A.y, cit->B.y, cit->C.y});
        std::tie(min.z, max.z) = std::minmax({min.z, max.z, cit->A.z, cit->B.z, cit->C.z});
    }
    min.x -= 1; min.y -= 1; min.z -= 1;
    max.x += 1; max.y += 1; max.z += 1;
    const point3 edge(max.x-min.x, max.y-min.y, max.z-min.z);
    const double delta = std::pow(edge.x*edge.y*edge.z/cell_count, 1.0/3);
    const index3 dim(std::ceil(edge.x/delta), std::ceil(edge.y/delta), std::ceil(edge.z/delta));
    trigrid triangle_grid(min, max, dim);
    triangle_grid.push(triangle_mesh);
    gstruct.org = make_float3(min.x, min.y, min.z);
    gstruct.dim = make_int4(dim.x, dim.y, dim.z, 0);
    gstruct.cell = make_float3(edge.x/dim.x, edge.y/dim.y, edge.z/dim.z);
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&gstruct.map_dev_p, gstruct.dim.x*gstruct.dim.y*gstruct.dim.z*sizeof(int));
        cuda_mem_scope<int>(gstruct.map_dev_p, gstruct.dim.x*gstruct.dim.y*gstruct.dim.z, [&](int* map_p) {
            for(int i = 0; i < gstruct.dim.x*gstruct.dim.y*gstruct.dim.z; i++)
                map_p[i] = -1;
        });
    });
    int2 size = make_int2(0, 0);
    triangle_grid.for_each_cell([&](const index3& i, const trimesh& cell_mesh) {
        size.x = std::max(size.x, cell_mesh.size());
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
                            triangle_grid.for_each_cell([&](const index3& i, const trimesh& cell_mesh) {
                                for(int x = 0; x < cell_mesh.size(); x++) {
                                    in_p[y][x] = cell_mesh[x].in;
                                    out_p[y][x] = cell_mesh[x].out;
                                    Ax_p[y][x] = cell_mesh[x].A.x;
                                    Ay_p[y][x] = cell_mesh[x].A.y;
                                    Az_p[y][x] = cell_mesh[x].A.z;
                                    Bx_p[y][x] = cell_mesh[x].B.x;
                                    By_p[y][x] = cell_mesh[x].B.y;
                                    Bz_p[y][x] = cell_mesh[x].B.z;
                                    Cx_p[y][x] = cell_mesh[x].C.x;
                                    Cy_p[y][x] = cell_mesh[x].C.y;
                                    Cz_p[y][x] = cell_mesh[x].C.z;
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