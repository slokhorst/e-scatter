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
#include <stack>
#include <vector>
#include <common/cuda_make_ptr.cuh>
#include <common/cuda_mem_scope.cuh>
#include <common/cuda_safe_call.cuh>

#include <iostream>

__host__ cuda_geometry_struct cuda_geometry_struct::create(const octree& root) {
    cuda_geometry_struct gstruct;
    gstruct.octree_dev_p = nullptr;
    gstruct.material_idx_in_dev_p = gstruct.material_idx_out_dev_p = nullptr;
    gstruct.triangle_r0x_dev_p = gstruct.triangle_r0y_dev_p = gstruct.triangle_r0z_dev_p = nullptr;
    gstruct.triangle_e1x_dev_p = gstruct.triangle_e1y_dev_p = gstruct.triangle_e1z_dev_p = nullptr;
    gstruct.triangle_e2x_dev_p = gstruct.triangle_e2y_dev_p = gstruct.triangle_e2z_dev_p = nullptr;
    if(root.triangles().empty())
        return gstruct;

    // sort octree nodes by location code (morton order)
    std::map<uint64_t,const octree*> morton_map;
    std::stack<const octree*> node_p_stack;
    node_p_stack.push(&root);
    while(!node_p_stack.empty()) {
        const octree* node_p = node_p_stack.top();
        node_p_stack.pop();
        morton_map[node_p->location()] = node_p;
        for(int octant = 0; octant < 8; octant++) {
            const octree* child_p = node_p->traverse(octant);
            if(child_p != nullptr)
                node_p_stack.push(child_p);
        }
    }

    // map triangles from octree to indices following location code
    std::vector<const triangle*> triangle_p_vec;
    std::map<const triangle*,int> triangle_p_map;
    for(auto morton_cit = morton_map.cbegin(); morton_cit != morton_map.cend(); morton_cit++) {
        const octree* node_p = morton_cit->second;
        if(node_p->is_leaf())
            for(const triangle* triangle_p : node_p->triangles())
                if(triangle_p_map.count(triangle_p) == 0) {
                    const int index = triangle_p_vec.size();
                    triangle_p_map[triangle_p] = index;
                    triangle_p_vec.push_back(triangle_p);
                }
    }

    // build linearized octree index table
    //  index=0 : child does not exist
    //  index>0 : non-leaf child with node indices
    //  index<0 : leaf child with triangle indices (index -1 means no triangle)
    std::vector<int> octree_vec;
    std::map<const octree*,int> node_p_map;
    for(auto morton_cit = morton_map.cbegin(); morton_cit != morton_map.cend(); morton_cit++) {
        const octree* node_p = morton_cit->second;
        const int index = octree_vec.size();
        node_p_map[node_p] = index;
        if(node_p->is_leaf()) {
            for(const triangle* triangle_p : node_p->triangles())
                octree_vec.push_back(triangle_p_map[triangle_p]);
            octree_vec.push_back(-1);
        } else {
            for(int octant = 0; octant < 8; octant++)
                octree_vec.push_back(0);
        }
    }
    for(auto cit = node_p_map.cbegin(); cit != node_p_map.cend(); cit++) {
        const octree* node_p = cit->first;
        const int index = cit->second;
        if(!node_p->is_leaf())
            for(int octant = 0; octant < 8; octant++) {
                const octree* child_p = node_p->traverse(octant);
                if(child_p != nullptr) {
                    octree_vec[index+octant] = node_p_map[child_p];
                    if(child_p->is_leaf())
                        octree_vec[index+octant] *= -1;
                }
            }
    }

    // copy octree to device memory
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&gstruct.octree_dev_p, octree_vec.size()*sizeof(int));
        cudaMemcpy(gstruct.octree_dev_p, octree_vec.data(), octree_vec.size()*sizeof(int), cudaMemcpyHostToDevice);
        gstruct.AABB_center = make_float3(root.center().x, root.center().y, root.center().z);
        gstruct.AABB_halfsize = make_float3(root.halfsize().x, root.halfsize().y, root.halfsize().z);
        gstruct.occupancy = root.occupancy();
    });

    // copy triangles to device memory
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&gstruct.material_idx_in_dev_p, triangle_p_vec.size()*sizeof(int));
        cudaMalloc(&gstruct.material_idx_out_dev_p, triangle_p_vec.size()*sizeof(int));
        cuda_mem_scope<int>(gstruct.material_idx_in_dev_p, triangle_p_vec.size(), [&](int* in_p) {
        cuda_mem_scope<int>(gstruct.material_idx_out_dev_p, triangle_p_vec.size(), [&](int* out_p) {
            for(size_t i = 0; i < triangle_p_vec.size(); i++) {
                in_p[i] = triangle_p_vec[i]->in;
                out_p[i] = triangle_p_vec[i]->out;
            }
        });
        });
    });
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&gstruct.triangle_r0x_dev_p, triangle_p_vec.size()*sizeof(float));
        cudaMalloc(&gstruct.triangle_r0y_dev_p, triangle_p_vec.size()*sizeof(float));
        cudaMalloc(&gstruct.triangle_r0z_dev_p, triangle_p_vec.size()*sizeof(float));
        cuda_mem_scope<float>(gstruct.triangle_r0x_dev_p, triangle_p_vec.size(), [&](float* r0x_p) {
        cuda_mem_scope<float>(gstruct.triangle_r0y_dev_p, triangle_p_vec.size(), [&](float* r0y_p) {
        cuda_mem_scope<float>(gstruct.triangle_r0z_dev_p, triangle_p_vec.size(), [&](float* r0z_p) {
            for(size_t i = 0; i < triangle_p_vec.size(); i++) {
                r0x_p[i] = triangle_p_vec[i]->A.x;
                r0y_p[i] = triangle_p_vec[i]->A.y;
                r0z_p[i] = triangle_p_vec[i]->A.z;
            }
        });
        });
        });
    });
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&gstruct.triangle_e1x_dev_p, triangle_p_vec.size()*sizeof(float));
        cudaMalloc(&gstruct.triangle_e1y_dev_p, triangle_p_vec.size()*sizeof(float));
        cudaMalloc(&gstruct.triangle_e1z_dev_p, triangle_p_vec.size()*sizeof(float));
        cuda_mem_scope<float>(gstruct.triangle_e1x_dev_p, triangle_p_vec.size(), [&](float* e1x_p) {
        cuda_mem_scope<float>(gstruct.triangle_e1y_dev_p, triangle_p_vec.size(), [&](float* e1y_p) {
        cuda_mem_scope<float>(gstruct.triangle_e1z_dev_p, triangle_p_vec.size(), [&](float* e1z_p) {
            for(size_t i = 0; i < triangle_p_vec.size(); i++) {
                e1x_p[i] = triangle_p_vec[i]->B.x-triangle_p_vec[i]->A.x;
                e1y_p[i] = triangle_p_vec[i]->B.y-triangle_p_vec[i]->A.y;
                e1z_p[i] = triangle_p_vec[i]->B.z-triangle_p_vec[i]->A.z;
            }
        });
        });
        });
    });
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&gstruct.triangle_e2x_dev_p, triangle_p_vec.size()*sizeof(float));
        cudaMalloc(&gstruct.triangle_e2y_dev_p, triangle_p_vec.size()*sizeof(float));
        cudaMalloc(&gstruct.triangle_e2z_dev_p, triangle_p_vec.size()*sizeof(float));
        cuda_mem_scope<float>(gstruct.triangle_e2x_dev_p, triangle_p_vec.size(), [&](float* e2x_p) {
        cuda_mem_scope<float>(gstruct.triangle_e2y_dev_p, triangle_p_vec.size(), [&](float* e2y_p) {
        cuda_mem_scope<float>(gstruct.triangle_e2z_dev_p, triangle_p_vec.size(), [&](float* e2z_p) {
            for(size_t i = 0; i < triangle_p_vec.size(); i++) {
                e2x_p[i] = triangle_p_vec[i]->C.x-triangle_p_vec[i]->A.x;
                e2y_p[i] = triangle_p_vec[i]->C.y-triangle_p_vec[i]->A.y;
                e2z_p[i] = triangle_p_vec[i]->C.z-triangle_p_vec[i]->A.z;
            }
        });
        });
        });
    });

    return gstruct;
}

__host__ void cuda_geometry_struct::release(cuda_geometry_struct& gstruct) {
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaFree(gstruct.octree_dev_p);
    });
    gstruct.octree_dev_p = nullptr;
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        for(int* p : {gstruct.material_idx_in_dev_p, gstruct.material_idx_out_dev_p})
            cudaFree(p);
        for(float* p : {gstruct.triangle_r0x_dev_p, gstruct.triangle_r0y_dev_p, gstruct.triangle_r0z_dev_p})
            cudaFree(p);
        for(float* p : {gstruct.triangle_e1x_dev_p, gstruct.triangle_e1y_dev_p, gstruct.triangle_e1z_dev_p})
            cudaFree(p);
        for(float* p : {gstruct.triangle_e2x_dev_p, gstruct.triangle_e2y_dev_p, gstruct.triangle_e2z_dev_p})
            cudaFree(p);
    });
    gstruct.material_idx_in_dev_p = gstruct.material_idx_out_dev_p = nullptr;
    gstruct.triangle_r0x_dev_p = gstruct.triangle_r0y_dev_p = gstruct.triangle_r0z_dev_p = nullptr;
    gstruct.triangle_e1x_dev_p = gstruct.triangle_e1y_dev_p = gstruct.triangle_e1z_dev_p = nullptr;
    gstruct.triangle_e2x_dev_p = gstruct.triangle_e2y_dev_p = gstruct.triangle_e2z_dev_p = nullptr;
}