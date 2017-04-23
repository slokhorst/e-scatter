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
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

#include <boost/program_options.hpp>
#include <curand_kernel.h>
#include <cub/cub.cuh>

#include "../common/archive.hh"
#include "../common/constant.hh"
#include "../common/profile_scope.hh"
#include "../cuda_common/cuda_mem_scope.cuh"
#include "../cuda_common/cuda_safe_call.cuh"
#include "../common/material.hh"

#include "cuda_kernels.cuh"
#include "octree.hh"

namespace po = boost::program_options;

uint64_t make_morton(const uint16_t x, const uint16_t y, const uint16_t z) {
    uint64_t morton = 0;
    for(int i = 0; i < 16; i++) {
        morton |= (x&(static_cast<uint64_t>(1)<<i))<<(2*i+0);
        morton |= (y&(static_cast<uint64_t>(1)<<i))<<(2*i+1);
        morton |= (z&(static_cast<uint64_t>(1)<<i))<<(2*i+2);
    }
    return morton;
}

std::unique_ptr<octree> load_octree_from_file(const std::string& file) {
    std::vector<triangle> triangle_vec;

    std::ifstream ifs(file);
    if(!ifs.is_open())
        return nullptr;
    int line_num = 1;
    while(!ifs.eof()) {
        std::string line_str;
        std::getline(ifs, line_str);
        line_num++;
        if(line_str.empty())
            continue;

        const size_t i = line_str.find_first_not_of(" \t");
        if(i < line_str.size())
        if(line_str[i] == '#')
            continue;

        std::stringstream ss;
        ss << line_str;
        std::vector<std::string> tag_vec;
        while(!ss.eof()) {
            std::string tag_str;
            ss >> tag_str;
            if(!tag_str.empty())
                tag_vec.push_back(tag_str);
        }
        if(tag_vec.size() != 11) {
            std::ostringstream oss;
            oss << "invalid number of columns in line '" << line_num << "'";
            throw std::runtime_error(oss.str());
        }
        int in, out;
        in = std::stoi(tag_vec[0]);
        out = std::stoi(tag_vec[1]);
        point3 A, B, C;
        A.x = std::stod(tag_vec[2]); A.y = std::stod(tag_vec[3]); A.z = std::stod(tag_vec[4]);
        B.x = std::stod(tag_vec[5]); B.y = std::stod(tag_vec[6]); B.z = std::stod(tag_vec[7]);
        C.x = std::stod(tag_vec[8]); C.y = std::stod(tag_vec[9]); C.z = std::stod(tag_vec[10]);
        triangle_vec.push_back(triangle(A, B, C, in, out));
    }
    ifs.close();

    if(triangle_vec.empty())
        return nullptr;

    point3 AABB_min, AABB_max;
    auto cit = triangle_vec.cbegin();
    std::tie(AABB_min.x, AABB_max.x) = std::minmax({cit->A.x, cit->B.x, cit->C.x});
    std::tie(AABB_min.y, AABB_max.y) = std::minmax({cit->A.y, cit->B.y, cit->C.y});
    std::tie(AABB_min.z, AABB_max.z) = std::minmax({cit->A.z, cit->B.z, cit->C.z});
    while(++cit != triangle_vec.cend()) {
        std::tie(AABB_min.x, AABB_max.x) =
            std::minmax({AABB_min.x, AABB_max.x, cit->A.x, cit->B.x, cit->C.x});
        std::tie(AABB_min.y, AABB_max.y) =
            std::minmax({AABB_min.y, AABB_max.y, cit->A.y, cit->B.y, cit->C.y});
        std::tie(AABB_min.z, AABB_max.z) =
            std::minmax({AABB_min.z, AABB_max.z, cit->A.z, cit->B.z, cit->C.z});
    }
    AABB_min -= point3(1, 1, 1);
    AABB_max += point3(1, 1, 1);

    std::unique_ptr<octree> octree_p(new octree(AABB_min, AABB_max));
    for(auto cit = triangle_vec.cbegin(); cit != triangle_vec.cend(); cit++)
        octree_p->insert(*cit);
    return octree_p;
}

int main(const int argc, char* argv[]) {
    double batch_factor;
    bool dry_run_flag;
    int prescan_size, cu_seed;
    std::string geometry_file;
    std::string particle_file;
    std::vector<std::string> material_file_vec;
    scatter_options opt;
    opt.generate_secondary_flag = true;

    std::string usage = "Usage: cdsem <geometry.tri> <particles.pri> [material0.mat] [...] [materialN.mat] [options]";
    po::options_description visible_options("available options");
    visible_options.add_options()
        ("help,h", "produce help message")
        ("batch-factor", po::value<double>(&batch_factor)->default_value(0.9),
            "fraction of the maximum possible batch size to use")
        ("dry-run", po::value<bool>(&dry_run_flag)->default_value(false),
            "no output will be generated")
        /*("random-seed", po:value<int>(&cu_seed)->default_value(-1),
            "give a random seed for cuda, -1 (default) uses timer.")*/
        ("prescan-size", po::value<int>(&prescan_size)->default_value(1000),
            "number of primary electrons used in the pre-scan")
        ("acoustic-phonon-loss", po::value<bool>(&opt.acoustic_phonon_loss_flag)->default_value(true),
            "enable acoustic phonon loss")
        ("optical-phonon-loss", po::value<bool>(&opt.optical_phonon_loss_flag)->default_value(true),
            "enable optical phonon loss")
        ("atomic-recoil-loss", po::value<bool>(&opt.atomic_recoil_loss_flag)->default_value(true),
            "enable atomic recoil loss")
        ("generate-secondary", po::value<bool>(&opt.generate_secondary_flag)->default_value(true),
            "enable or disable creation of secondary electrons")
        ("instantaneous-momentum", po::value<bool>(&opt.instantaneous_momentum_flag)->default_value(true),
            "use a random instantaneous momentum of the secondary prior to collision?")
        ("momentum-conservation", po::value<bool>(&opt.momentum_conservation_flag)->default_value(true),
            "enable momentum conservation. note: without momentum conservation -> forward scattering of the incident electron is assumed")
        ("quantum-transmission", po::value<bool>(&opt.quantum_transmission_flag)->default_value(true),
            "enable quantum transmission")
        ("interface-refraction", po::value<bool>(&opt.interface_refraction_flag)->default_value(true),
            "enable interface refraction")
        ("interface-absorption", po::value<bool>(&opt.interface_absorption_flag)->default_value(false),
            "enable interface absorption for compliance with Kieft")
    ;
    po::options_description hidden_options;
    hidden_options.add_options()
        ("geometry-file", po::value<std::string>(&geometry_file))
        ("particle-file", po::value<std::string>(&particle_file))
        ("material-file", po::value<std::vector<std::string>>(&material_file_vec))
    ;
    po::positional_options_description posoptions;
    posoptions.add("geometry-file", 1);
    posoptions.add("particle-file", 1);
    posoptions.add("material-file",-1);
    po::options_description options;
    options.add(visible_options).add(hidden_options);

    try {
        po::variables_map vars;
        po::store(po::command_line_parser(argc, argv).options(options).positional(posoptions).run(), vars);
        if(vars.count("help")) {
            std::clog << usage << std::endl
                      << visible_options << std::endl;
            return EXIT_SUCCESS;
        }
        po::notify(vars);

        if(vars.count("geometry-file")<1)
            throw std::runtime_error("no geometry file defined");
        if(vars.count("particle-file")<1)
            throw std::runtime_error("no particle file defined");

    } catch(const std::exception& e) {
        std::clog << e.what() << std::endl;
        std::clog << usage << std::endl
                  << visible_options << std::endl;
        return EXIT_FAILURE;
    }

    std::clog << " >> loading geometry file='" << geometry_file << "'";
    std::clog << std::flush;
    std::unique_ptr<octree> octree_p = load_octree_from_file(geometry_file);
    if(octree_p == nullptr) {
        std::clog << std::endl;
        throw std::runtime_error("no geometry");
    }
    std::clog << " triangles=" << octree_p->triangles().size();
    std::clog << " tree_depth=" << octree_p->depth();
    std::clog << " occupancy=" << octree_p->occupancy();
    std::clog << " min=(" << octree_p->center().x-octree_p->halfsize().x+1;
    std::clog << "," << octree_p->center().y-octree_p->halfsize().y+1;
    std::clog << "," << octree_p->center().z-octree_p->halfsize().z+1 << ")";
    std::clog << " max=(" << octree_p->center().x+octree_p->halfsize().x-1;
    std::clog << "," << octree_p->center().y+octree_p->halfsize().y-1;
    std::clog << "," << octree_p->center().z+octree_p->halfsize().z-1 << ")";
    std::map<int,int> material_map;
    for(const triangle* triangle_p : octree_p->triangles()) {
        material_map[triangle_p->in]++;
        material_map[triangle_p->out]++;
    }
    for(auto cit = material_map.cbegin(); cit != material_map.cend(); cit++) {
        std::clog << " idx:cnt=" << cit->first << ":" << cit->second;
        if(cit->first > (int)material_file_vec.size()-1)
            throw std::runtime_error("no material file specified with idx="+std::to_string(cit->first));
    }
    std::clog << std::endl;

    std::clog << " >> loading particles file='" << particle_file << "'";
    std::clog << std::flush;
    struct particle {
        float3 pos;
        float3 dir;
        float K;
        int tag;
    };
    std::vector<particle> particle_vec;
    std::vector<std::pair<int,int>> tag_map;
    std::ifstream ifs(particle_file, std::ifstream::binary);
    if(!ifs.is_open()) {
        std::clog << std::endl;
        throw std::ios_base::failure("failed to open '"+particle_file+"' for reading");
    }
    while(!ifs.eof()) {
        particle primary;
        ifs.read(reinterpret_cast<char*>(&(primary.pos.x)), sizeof(primary.pos.x));
        ifs.read(reinterpret_cast<char*>(&(primary.pos.y)), sizeof(primary.pos.y));
        ifs.read(reinterpret_cast<char*>(&(primary.pos.z)), sizeof(primary.pos.z));
        ifs.read(reinterpret_cast<char*>(&(primary.dir.x)), sizeof(primary.dir.x));
        ifs.read(reinterpret_cast<char*>(&(primary.dir.y)), sizeof(primary.dir.y));
        ifs.read(reinterpret_cast<char*>(&(primary.dir.z)), sizeof(primary.dir.z));
        ifs.read(reinterpret_cast<char*>(&primary.K), sizeof(primary.K));
        int2 pixel;
        ifs.read(reinterpret_cast<char*>(&(pixel.x)), sizeof(pixel.x));
        ifs.read(reinterpret_cast<char*>(&(pixel.y)), sizeof(pixel.y));
        primary.tag = tag_map.size();
        if(ifs.eof())
            break;
        particle_vec.push_back(primary);
        tag_map.push_back(std::make_pair(pixel.x, pixel.y));
    }
    ifs.close();
    std::clog << " count=" << particle_vec.size();
    std::clog << std::endl;
    const int particle_cnt = particle_vec.size();
    if(particle_cnt == 0)
        throw std::runtime_error("no particles");

    std::clog << " >> sorting particles";
    std::clog << std::flush;
    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(particle_vec.begin(), particle_vec.end(), std::default_random_engine(seed));
    if(prescan_size < particle_cnt) {
        std::sort(particle_vec.begin()+prescan_size, particle_vec.end(), [&](const particle& p1, const particle& p2) {
            const point3 AABB_min = octree_p->center()-octree_p->halfsize();
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
    }
    std::clog << std::endl;

    std::vector<material> material_vec;
    for(const std::string& material_file : material_file_vec) {
        std::clog << " >> loading material";
        std::clog << " index=" << material_vec.size();
        std::clog << " file='" << material_file << "'";
        std::clog << std::flush;
        std::ifstream ifs(material_file, std::ifstream::binary);
        if(!ifs.is_open()) {
            std::clog << std::endl;
            throw std::ios_base::failure("failed to open '"+material_file+"' for reading");
        }
        archive::istream ia(ifs);
        material_vec.push_back(material());
        ia >> material_vec.back();
        ifs.close();
        std::clog << " fermi=" << material_vec.back().fermi()/constant::ec;
        std::clog << " barrier=" << material_vec.back().barrier()/constant::ec;
        if(material_vec.back().band_gap().is_defined())
            std::clog << " band_gap=" << material_vec.back().band_gap()()/constant::ec;
        std::clog << " phonon_loss=" << material_vec.back().phonon_loss()/constant::ec;
        std::clog << std::endl;
    }
    if(material_map.crbegin()->first > static_cast<int>(material_vec.size()))
        throw std::runtime_error("incomplete material set");


    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    });


    const int max_warp_count = 31250/5;
    const int warps_per_block = 4;
    const int threads_per_warp = 32;
    const int threads_per_block = warps_per_block*threads_per_warp;

    const int capacity = max_warp_count*threads_per_warp;

    std::clog << " >> pushing geometry to device";
    std::clog << std::flush;
    cuda_geometry_struct gstruct(cuda_geometry_struct::create(*octree_p));
    std::clog << std::endl;

    cuda_particle_struct pstruct(cuda_particle_struct::create(capacity));
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

    cuda_material_struct mstruct(cuda_material_struct::create(material_vec.size()));
    for(int i = 0; i < mstruct.capacity; i++)
        mstruct.assign(i, material_vec[i]);

    int* radix_index_dev_p = nullptr;
    void* radix_dump_dev_p = nullptr;
    void* radix_temp_dev_p = nullptr;
    size_t radix_temp_size = 0;
    curandState* rand_state_dev_p = nullptr;

    /* seed = (cu_seed == -1 ?
               std::chrono::system_clock::now().time_since_epoch().count()
             : cu_seed); */

    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&radix_index_dev_p, capacity*sizeof(int));
        cuda_mem_scope<int>(radix_index_dev_p, capacity, [&](int* index_p) {
            for(int i = 0; i < capacity; i++)
                index_p[i] = i;
        });
        cudaMalloc(&radix_dump_dev_p, capacity*sizeof(uint8_t));
        cub::DeviceRadixSort::SortPairs<uint8_t,int>(
            nullptr, radix_temp_size,
            pstruct.status_dev_p, static_cast<uint8_t*>(radix_dump_dev_p),
            radix_index_dev_p, pstruct.particle_idx_dev_p,
            capacity
        );
        radix_temp_size++;
        cudaMalloc(&radix_temp_dev_p, radix_temp_size);
        cudaMalloc(&rand_state_dev_p, capacity*sizeof(curandState));
        std::clog << " >> initializing random states";
        std::clog << std::flush;
        cuda_init_rand_state<<<1+capacity/threads_per_block,threads_per_block>>>(
            rand_state_dev_p, seed, capacity, opt
        );
        cudaDeviceSynchronize();
        std::clog << std::endl;
    });

    auto __execute_iteration = [&] {
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_init_trajectory<<<1+capacity/threads_per_block,threads_per_block>>>(
                pstruct, gstruct, mstruct, rand_state_dev_p, opt
            );
            cuda_update_trajectory<<<1+capacity/threads_per_block,threads_per_block>>>(
                pstruct, gstruct, mstruct, opt
            );
            cub::DeviceRadixSort::SortPairs<uint8_t,int>(
                radix_temp_dev_p, radix_temp_size,
                pstruct.status_dev_p, static_cast<uint8_t*>(radix_dump_dev_p),
                radix_index_dev_p, pstruct.particle_idx_dev_p,
                capacity, 0, 2
            );
            cuda_intersection_event<<<1+capacity/threads_per_block,threads_per_block>>>(
                pstruct, gstruct, mstruct, rand_state_dev_p, opt
            );
            cuda_elastic_event<<<1+capacity/threads_per_block,threads_per_block>>>(
                pstruct, mstruct, rand_state_dev_p, opt
            );
            cuda_inelastic_event<<<1+capacity/threads_per_block,threads_per_block>>>(
                pstruct, mstruct, rand_state_dev_p, opt
            );
        });
    };

    auto __flush_detected_particles = [&](std::ostream& os) {
        pstruct.for_each(cuda_particle_struct::DETECTED, [&](float3 pos, float3 dir, float K, int tag) {
            if(!dry_run_flag) {
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
            }
        });
        pstruct.flush();
    };

    int batch_index = pstruct.push(pos_vec.data(), dir_vec.data(), K_vec.data(), tag_vec.data(), prescan_size);
    std::vector<std::pair<int,int>> prescan_stats_vec;
    prescan_stats_vec.push_back(std::make_pair(batch_index, 0));
    while(prescan_stats_vec.back().first > 0) {
        __execute_iteration();
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_mem_scope<uint8_t>(pstruct.status_dev_p, capacity, [&](uint8_t* status_p) {
                int running_count = 0;
                int detected_count = 0;
                for(int i = 0; i < capacity; i++)
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
        std::clog << " \r";
        std::clog << " >> executing pre-scan";
        std::clog << " running_count=" << prescan_stats_vec.back().first;
        std::clog << " detected_count=" << prescan_stats_vec.back().second;
        std::clog << std::flush;
    }
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
    const int batch_size = batch_factor*capacity/accumulator;

    const size_t time_lapse = profile_scope([&]{
        int running_count;
        do {
            const int push_count = std::min(particle_cnt-batch_index, batch_size);
            if(push_count > 0)
                batch_index += pstruct.push(pos_vec.data()+batch_index, dir_vec.data()+batch_index, K_vec.data()+batch_index, tag_vec.data()+batch_index, push_count);
            for(int i = 0; i < frame_size; i++) {
                __execute_iteration();
                cudaDeviceSynchronize();
            }
            running_count = 0;
            cuda_safe_call(__FILE__, __LINE__, [&]() {
                cuda_mem_scope<uint8_t>(pstruct.status_dev_p, capacity, [&](uint8_t* status_p) {
                    for(int i = 0; i < capacity; i++)
                        switch(status_p[i]) {
                            case cuda_particle_struct::TERMINATED:
                            break;
                            default:
                                running_count++;
                            break;
                        }
                }); // status_dev_p
            });
            std::clog << " \r";
            std::clog << " >> executing exposure";
            std::clog << " pct=" << 100*batch_index/(particle_cnt-1) << "%";
            std::clog << " frame_size=" << frame_size;
            std::clog << " batch_size=" << batch_size;
            std::clog << " running_count=" << running_count;
            std::clog << std::flush;
            __flush_detected_particles(std::cout);
        } while(running_count > 0);
        std::clog << std::endl;
    });
    std::clog << " >> time_lapse=" << time_lapse << std::endl;

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
