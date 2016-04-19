/**
 * @file src/cdsem/main.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <algorithm>
#include <cinttypes>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <curand_kernel.h>
#include <cub/cub.cuh>
#include <common/archive.hh>
#include <common/profile_scope.hh>
#include <common/cuda_mem_scope.cuh>
#include <common/cuda_safe_call.cuh>
#include "cuda_kernels.cuh"
#include "material.hh"
#include "octree.hh"

class sc {
public:
	sc() : _cursor_vec(std::vector<char>{ '|', '/', '-', '\\' }) {}
    void reset() {
       _cursor_index = 0;
    }
    void print() {
        std::clog << "    \r [" << _cursor_vec[(_cursor_index++)%_cursor_vec.size()] << "]";
    }
    void finish() const {
        std::clog << "    \r [*]";
    }
private:
    std::vector<char> _cursor_vec;
    size_t _cursor_index = 0;
} spin_cursor;

uint64_t make_morton(const uint16_t x, const uint16_t y, const uint16_t z) {
    uint64_t morton = 0;
    for(int i = 0; i < 16; i++) {
        morton |= (x&(static_cast<uint64_t>(1)<<i))<<(2*i+0);
        morton |= (y&(static_cast<uint64_t>(1)<<i))<<(2*i+1);
        morton |= (z&(static_cast<uint64_t>(1)<<i))<<(2*i+2);
    }
    return morton;
}

int main(const int argc, char* argv[]) {
    if(argc < 3)
        return EXIT_FAILURE;

    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    });

    const int electron_capacity = 1e6;
    const int prescan_size = 256;

    std::vector<triangle> triangle_mesh;
    std::ifstream tri_ifs(argv[1]);
    while(!tri_ifs.eof()) {
        int in, out;
        tri_ifs >> in >> out;
        if((triangle_mesh.size()%10240 == 0) || (tri_ifs.eof())) {
            spin_cursor.print();
            std::clog << " loading triangles";
            std::clog << " file='" << argv[1] << "'";
            std::clog << " count=" << triangle_mesh.size();
            std::clog << std::flush;
        }
        if(tri_ifs.eof())
            break;
        point3 A, B, C;
        tri_ifs >> A.x >> A.y >> A.z;
        tri_ifs >> B.x >> B.y >> B.z;
        tri_ifs >> C.x >> C.y >> C.z;
        triangle_mesh.push_back(triangle(A, B, C, in, out));
    }
    spin_cursor.finish();
    std::clog << std::endl;
    tri_ifs.close();

    point3 AABB_min, AABB_max;
    AABB_min = AABB_max = triangle_mesh.front().A;
    for(auto cit = triangle_mesh.cbegin(); cit != triangle_mesh.cend(); cit++) {
        AABB_min.x = std::min({AABB_min.x, cit->A.x, cit->B.x, cit->C.x});
        AABB_min.y = std::min({AABB_min.y, cit->A.y, cit->B.y, cit->C.y});
        AABB_min.z = std::min({AABB_min.z, cit->A.z, cit->B.z, cit->C.z});
        AABB_max.x = std::max({AABB_max.x, cit->A.x, cit->B.x, cit->C.x});
        AABB_max.y = std::max({AABB_max.y, cit->A.y, cit->B.y, cit->C.y});
        AABB_max.z = std::max({AABB_max.z, cit->A.z, cit->B.z, cit->C.z});
    }
    std::clog << " [*] axis aligned bounding box";
    std::clog << " min=(" << AABB_min.x << "," << AABB_min.y << "," << AABB_min.z << ")";
    std::clog << " max=(" << AABB_max.x << "," << AABB_max.y << "," << AABB_max.z << ")";
    std::clog << std::endl;

    octree triangle_tree(AABB_min-point3(1, 1, 1), AABB_max+point3(1, 1, 1));
    spin_cursor.reset();
    for(size_t i = 0; i < triangle_mesh.size(); i++) {
        if(i%10240 == 0) {
            spin_cursor.print();
            std::clog << " building octree";
        }
        triangle_tree.insert(triangle_mesh[i]);
    }
    std::clog << " triangle_count=" << triangle_tree.count();
    std::clog << " tree_depth=" << triangle_tree.depth();
    std::clog << " occupancy=" << triangle_tree.occupancy();
    spin_cursor.finish();
    std::clog << std::endl;

    std::clog << " [*] initializing octree based geometry";
    cuda_geometry_struct gstruct(cuda_geometry_struct::create(triangle_tree));
    std::clog << std::endl;

    std::vector<material> mat_vec;
    for(const std::string& mat_file : {"../data/silicon.mat", "../data/pmma.mat"}) {
        std::clog << " [*] loading material";
        std::clog << " index=" << mat_vec.size();
        std::clog << " file='" << mat_file << "'";
        std::clog << std::endl;
        std::ifstream ifs(mat_file, std::ifstream::in | std::ifstream::binary);
        if(!ifs.is_open())
            throw std::ios_base::failure("failed to open '"+mat_file+"' for reading");
        archive::istream ia(ifs);
        mat_vec.push_back(material());
        ia >> mat_vec.back();
        ifs.close();
    }

    cuda_material_struct mstruct(cuda_material_struct::create(mat_vec.size()));
    for(int i = 0; i < mstruct.capacity; i++)
        mstruct.assign(i, mat_vec[i]);

    cuda_particle_struct pstruct(cuda_particle_struct::create(electron_capacity));
    struct particle {
        float3 pos;
        float3 dir;
        float K;
        int tag;
    };
    std::vector<particle> particle_vec;
    std::map<int,std::pair<int,int>> tag_map;
    std::ifstream pri_ifs(argv[2]);
    spin_cursor.reset();
    while(!pri_ifs.eof()) {
        particle primary;
        pri_ifs.read(reinterpret_cast<char*>(&(primary.pos.x)), sizeof(primary.pos.x));
        pri_ifs.read(reinterpret_cast<char*>(&(primary.pos.y)), sizeof(primary.pos.y));
        pri_ifs.read(reinterpret_cast<char*>(&(primary.pos.z)), sizeof(primary.pos.z));
        pri_ifs.read(reinterpret_cast<char*>(&(primary.dir.x)), sizeof(primary.dir.x));
        pri_ifs.read(reinterpret_cast<char*>(&(primary.dir.y)), sizeof(primary.dir.y));
        pri_ifs.read(reinterpret_cast<char*>(&(primary.dir.z)), sizeof(primary.dir.z));
        pri_ifs.read(reinterpret_cast<char*>(&primary.K), sizeof(primary.K));
        int2 pixel;
        pri_ifs.read(reinterpret_cast<char*>(&(pixel.x)), sizeof(pixel.x));
        pri_ifs.read(reinterpret_cast<char*>(&(pixel.y)), sizeof(pixel.y));
        primary.tag = tag_map.size();
        if((primary.tag%65536 == 0) || (pri_ifs.eof())) {
            spin_cursor.print();
            std::clog << " loading particles";
            std::clog << " file='" << argv[2] << "'";
            std::clog << " count=" << particle_vec.size();
            std::clog << std::flush;
        }
        if(pri_ifs.eof())
            break;
        if(primary.pos.x < -32)
            primary.pos.x = -64-primary.pos.x;
        if(primary.pos.x > 32)
            primary.pos.x = 64-primary.pos.x;
        if(primary.pos.y < -750)
            primary.pos.y = -1500-primary.pos.y;
        if(primary.pos.y > 750)
            primary.pos.y = 1500-primary.pos.y;
        particle_vec.push_back(primary);
        tag_map[primary.tag] = std::make_pair(pixel.x, pixel.y);
    }
    spin_cursor.finish();
    std::clog << std::endl;
    pri_ifs.close();
    const int np = particle_vec.size();
    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(particle_vec.begin(), particle_vec.end(), std::default_random_engine(seed));
    if(prescan_size < np) {
        std::clog << " [*] sorting particles" << std::flush;
        std::sort(particle_vec.begin()+prescan_size, particle_vec.end(), [&](const particle& p1, const particle& p2) {
            const uint16_t x1 = (p1.pos.x-AABB_min.x);
            const uint16_t y1 = (p1.pos.y-AABB_min.y);
            const uint16_t z1 = (p1.pos.z-AABB_min.z);
            const uint16_t x2 = (p2.pos.x-AABB_min.x);
            const uint16_t y2 = (p2.pos.y-AABB_min.y);
            const uint16_t z2 = (p2.pos.z-AABB_min.z);
            if(make_morton(x1, y1, z1) < make_morton(x2, y2, z2))
                return true;
            return false;
        });
        std::clog << std::endl;
    }

    std::vector<float3> pos_vec;
    std::vector<float3> dir_vec;
    std::vector<float> K_vec;
    std::vector<int> tag_vec;
    for(size_t i = 0; i < particle_vec.size(); i++) {
        pos_vec.push_back(particle_vec[i].pos);
        dir_vec.push_back(particle_vec[i].dir);
        K_vec.push_back(particle_vec[i].K);
        tag_vec.push_back(particle_vec[i].tag);
    }

    int* radix_index_dev_p = nullptr;
    void* radix_dump_dev_p = nullptr;
    void* radix_temp_dev_p = nullptr;
    size_t radix_temp_size = 0;
    curandState* rand_state_dev_p = nullptr;
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&radix_index_dev_p, electron_capacity*sizeof(int));
        cuda_mem_scope<int>(radix_index_dev_p, electron_capacity, [&](int* index_p) {
            for(int i = 0; i < electron_capacity; i++)
                index_p[i] = i;
        });
        cudaMalloc(&radix_dump_dev_p, electron_capacity*sizeof(uint8_t));
        cub::DeviceRadixSort::SortPairs<uint8_t,int>(
            nullptr, radix_temp_size,
            pstruct.status_dev_p, static_cast<uint8_t*>(radix_dump_dev_p),
            radix_index_dev_p, pstruct.particle_idx_dev_p,
            electron_capacity
        );
        radix_temp_size++;
        cudaMalloc(&radix_temp_dev_p, radix_temp_size);
        cudaMalloc(&rand_state_dev_p, electron_capacity*sizeof(curandState));
        cuda_init_rand_state<<<1+electron_capacity/128,128>>>(
            rand_state_dev_p, 0, electron_capacity
        );
    });

    auto __execute_iteration = [&] {
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_init_trajectory<<<1+electron_capacity/128,128>>>(pstruct, gstruct, mstruct, rand_state_dev_p);
        });
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_update_trajectory<<<1+electron_capacity/128,128>>>(pstruct, gstruct, mstruct);
        });
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cub::DeviceRadixSort::SortPairs<uint8_t,int>(
                radix_temp_dev_p, radix_temp_size,
                pstruct.status_dev_p, static_cast<uint8_t*>(radix_dump_dev_p),
                radix_index_dev_p, pstruct.particle_idx_dev_p,
                electron_capacity, 0, 2
            );
        });
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_apply_intersection_event<<<1+electron_capacity/128,128>>>(pstruct, gstruct, mstruct, rand_state_dev_p);
        });
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_apply_elastic_event<<<1+electron_capacity/128,128>>>(pstruct, mstruct, rand_state_dev_p);
        });
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_apply_inelastic_event<<<1+electron_capacity/128,128>>>(pstruct, mstruct, rand_state_dev_p);
        });
    };

    auto __flush_detected_particles = [&](std::ostream& os) {
        pstruct.for_each(cuda_particle_struct::DETECTED, [&](float3 pos, float3 dir, float K, int tag) {
            const int2 pixel = make_int2(tag_map[tag].first, tag_map[tag].second);
            os.write(reinterpret_cast<const char*>(&(pos.x)), sizeof(pos.x));
            os.write(reinterpret_cast<const char*>(&(pos.y)), sizeof(pos.y));
            os.write(reinterpret_cast<const char*>(&(pos.z)), sizeof(pos.z));
            os.write(reinterpret_cast<const char*>(&(dir.x)), sizeof(dir.x));
            os.write(reinterpret_cast<const char*>(&(dir.y)), sizeof(dir.y));
            os.write(reinterpret_cast<const char*>(&(dir.z)), sizeof(dir.z));
            os.write(reinterpret_cast<const char*>(&K), sizeof(K));
            os.write(reinterpret_cast<const char*>(&(pixel.x)), sizeof(pixel.x));
            os.write(reinterpret_cast<const char*>(&(pixel.y)), sizeof(pixel.y));
        });
        pstruct.flush();
    };

    int batch_index = pstruct.push(pos_vec.data(), dir_vec.data(), K_vec.data(), tag_vec.data(), prescan_size);
    std::vector<std::pair<int,int>> prescan_stats_vec;
    prescan_stats_vec.push_back(std::make_pair(batch_index, 0));
    spin_cursor.reset();
    while(prescan_stats_vec.back().first > 0) {
        __execute_iteration();
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_mem_scope<uint8_t>(pstruct.status_dev_p, electron_capacity, [&](uint8_t* status_p) {
                int running_count = 0;
                int detected_count = 0;
                for(int i = 0; i < electron_capacity; i++)
                    switch(status_p[i]) {
                        case cuda_particle_struct::DETECTED:
                            detected_count++;
                        break;
                        case cuda_particle_struct::TERMINATED:
                        break;
                        default:
                            running_count++;
                        break;
                    }
                prescan_stats_vec.push_back(std::make_pair(running_count, detected_count));
            }); // status_dev_p
        });
        spin_cursor.print();
        std::clog << " executing prescan";
        std::clog << " running_count=" << prescan_stats_vec.back().first;
        std::clog << " detected_count=" << prescan_stats_vec.back().second;
        std::clog << std::flush;
    }
    spin_cursor.finish();
    std::clog << std::endl;
    int frame_size = 1;
    for(size_t i = 1; i < prescan_stats_vec.size(); i++)
        if(prescan_stats_vec[i].first > prescan_stats_vec[frame_size].first)
            frame_size = i;
    float accumulator = 0;
    for(size_t i = 2*frame_size; i < prescan_stats_vec.size(); i += frame_size)
        accumulator += 1.0*prescan_stats_vec[i].first/prescan_size;
    accumulator += 2.0*prescan_stats_vec[frame_size].first/prescan_size;
    accumulator += 2.0*prescan_stats_vec[frame_size].second/prescan_size;
    int batch_size = 0.95*electron_capacity/accumulator;

    std::clog << " [*] time_lapse=" << profile_scope([&]{
        int running_count;
        spin_cursor.reset();
        do {
            batch_size = std::min(np-batch_index, batch_size);
            if(batch_size > 0)
                batch_index += pstruct.push(pos_vec.data()+batch_index, dir_vec.data()+batch_index, K_vec.data()+batch_index, tag_vec.data()+batch_index, batch_size);
            for(int i = 0; i < frame_size; i++)
                __execute_iteration();
            __flush_detected_particles(std::cout);
            running_count = 0;
            cuda_safe_call(__FILE__, __LINE__, [&]() {
                cuda_mem_scope<uint8_t>(pstruct.status_dev_p, electron_capacity, [&](uint8_t* status_p) {
                    for(int i = 0; i < electron_capacity; i++)
                        switch(status_p[i]) {
                            case cuda_particle_struct::TERMINATED:
                            break;
                            default:
                                running_count++;
                            break;
                        }
                }); // status_dev_p
            });
            spin_cursor.print();
            std::clog << " executing exposure";
            std::clog << " pct=" << 100*batch_index/(np-1) << "%";
            std::clog << " frame_size=" << frame_size;
            std::clog << " batch_size=" << batch_size;
            std::clog << " running_count=" << running_count;
            std::clog << std::flush;
        } while(running_count > 0);
        spin_cursor.finish();
        std::clog << std::endl;
    }) << std::endl;

    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaFree(radix_temp_dev_p);
        cudaFree(radix_dump_dev_p);
        cudaFree(radix_index_dev_p);
        cudaFree(rand_state_dev_p);
    });
    cuda_geometry_struct::release(gstruct);
    cuda_material_struct::release(mstruct);
    cuda_particle_struct::release(pstruct);

    return EXIT_SUCCESS;
}