/**
 * @file src/cdsem/cuda_kernels.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "cuda_kernels.cuh"
#include <cassert>
#include <cfloat>
#include <cuda_common/cuda_make_ptr.cuh>
#include <cuda_common/cuda_vec3_math.cuh>

__device__ const float _eps = 10.0f*FLT_EPSILON;
__device__ const float _me = 9.10938356e-31f; // electron mass [kg]
//__device__ const float _c2 = 8.98755179e+16f; // speed of light squared [m/s]
//__device__ const float _eV = 1.60217662e-19f; // one electron volt [C]
__device__ const float _NA = 6.02214086e+23f; // Avogadro's constant [/mol]
//__device__ const float _mc2 = _me*_c2/_eV;    // rest mass of electron [eV]
__device__ const float _pi = 3.14159265f;

__device__ float posf(float x) {
    return fmaxf(0.0f, x);
}

/*! \brief Scalar value clamping
 *  @param[in] x
 *  @param[in] x1
 *  @param[in] x2
 *  @return
 *   Clamped value such that x1<=x<=x2.
 */
template<typename T>
__device__ T clamp(T x, T x1, T x2) {
    return min(max(x, x1), x2);
}

/*! \brief Linear 1D interpolation routine.
 *  @param[in] ptr
 *   pointer to raw interpolation data table
 *  @param[in] pitch
 *   width in bytes
 *  @param[in] height
 *   number of rows
 *  @return
 *   Linearly interpolated value.
 */
template<typename T>
__device__ T interp1(const T* ptr, int pitch, int height, int3 offset, int dim, float x) {
    const int i = clamp(__float2int_rd(x*(dim-1)), 0, dim-2);
    const float s = x*(dim-1)-i;

    const T* p0 = cuda_make_ptr(ptr, pitch, height, offset.y, offset.z)+offset.x;
    return (1.0f-s)*p0[i]+s*p0[i+1];
}

/*! \brief Linear 2D interpolation routine.
 *  @param[in] ptr
 *   pointer to raw interpolation table
 *  @param[in] pitch
 *   width in bytes
 *  @param[in] height
 *   number of rows
 *  @return
 *   Linearly interpolated value.
 */
template<typename T>
__device__ T interp2(const T* ptr, int pitch, int height, int3 offset, int2 dim, float x, float y) {
    const int i = clamp(__float2int_rd(x*(dim.x-1)), 0, dim.x-2);
    const float s = x*(dim.x-1)-i;
    const int j = clamp(__float2int_rd(y*(dim.y-1)), 0, dim.y-2);
    const float t = y*(dim.y-1)-j;

    const T* p0 = cuda_make_ptr(ptr, pitch, height, offset.y+j+0, offset.z)+offset.x;
    const T* p1 = cuda_make_ptr(ptr, pitch, height, offset.y+j+1, offset.z)+offset.x;
    return (1.0f-s)*(1.0f-t)*p0[i]+s*(1.0f-t)*p0[i+1]+(1.0f-s)*t*p1[i]+s*t*p1[i+1];
}

__device__ float3 make_unit_vec(float cos_theta, float phi) {
    float sin_phi, cos_phi;
    sincosf(phi, &sin_phi, &cos_phi);
    cos_theta = clamp(cos_theta, -1.0f, 1.0f);
    const float sin_theta = sqrtf(1.0f-cos_theta*cos_theta);
    return make_float3(sin_theta*cos_phi, sin_theta*sin_phi, cos_theta);
}

/*! \brief Determination of a normal vector
 *  @return
 *   directional vector for which holds that dot_product(dir, make_normal_vec(dir, phi)) is zero.
 */
__device__ float3 make_normal_vec(float3 dir, float phi) {
    float sin_azimuth, cos_azimuth;
    sincosf(atan2f(dir.y, dir.x), &sin_azimuth, &cos_azimuth);

    const float3 unit_v = make_float3(
        dir.z*cos_azimuth,
        dir.z*sin_azimuth,
        -sqrtf(dir.x*dir.x+dir.y*dir.y)
    );
    const float3 unit_u = cross_product(unit_v, dir);

    float sin_phi, cos_phi;
    sincosf(phi, &sin_phi, &cos_phi);
    return unit_u*cos_phi+unit_v*sin_phi;
}

/*! \brief Inside AABB detection.
 *  @param[in] pos
 *  @param[in] center
 *  @param[in] halfsize
 *  @return
 */
__device__ bool inside_AABB(float3 pos, float3 center, float3 halfsize) {
    if((pos.x > center.x-halfsize.x) && (pos.x < center.x+halfsize.x))
    if((pos.y > center.y-halfsize.y) && (pos.y < center.y+halfsize.y))
    if((pos.z > center.z-halfsize.z) && (pos.z < center.z+halfsize.z))
        return true;
    return false;
}

__device__ float3 AABB_intersect(float3 pos, float3 dir, float3 center, float3 halfsize) {
    return make_float3(
        (center.x+copysignf(halfsize.x+_eps, dir.x)-pos.x)/dir.x,
        (center.y+copysignf(halfsize.y+_eps, dir.y)-pos.y)/dir.y,
        (center.z+copysignf(halfsize.z+_eps, dir.z)-pos.z)/dir.z
    );
}

__device__ float triangle_intersect(float3 pos, float3 dir, float3 O, float3 e1, float3 e2) {
    // T. Möller and B. Trumbore, Journal of Graphics Tools, 2(1):21--28, 1997.
    const float3 pvec = cross_product(dir, e2);
    const float det = dot_product(e1, pvec);
    if(fabsf(det) < _eps)
        return -1.0f;
    const float u = dot_product(pos-O, pvec)/det;
    if((u < -_eps) || (u > 1.0f+_eps))
        return -1.0f;
    const float3 qvec = cross_product(pos-O, e1);
    const float v = dot_product(dir, qvec)/det;
    if((v < -_eps) || (u+v > 1.0f+_eps))
        return -1.0f;
    const float t = dot_product(e2, qvec)/det;
    if(t < 0)
        return -1.0f;
    return t;
}

class curand_wrapper {
public:
    __device__ curand_wrapper(curandState& rand_state) {
        rand_state_dev_p = &rand_state;
        local_rand_state = rand_state;
    }
    __device__ curand_wrapper(const curand_wrapper&) = delete;
    __device__ ~curand_wrapper() {
        *rand_state_dev_p = local_rand_state;
    }
    __device__ curand_wrapper& operator=(const curand_wrapper&) = delete;
    __device__ float uniform() {
        return curand_uniform(&local_rand_state);
    }
private:
    curandState* rand_state_dev_p;
    curandState local_rand_state;
};

__global__ void cuda_init_rand_state(curandState* rand_state_p, unsigned long long seed, int n, scatter_options opt) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= n)
        return;
    curand_init(seed, i, 0, &(rand_state_p[i]));
}

__global__ void cuda_init_trajectory(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p, scatter_options opt) {
    const int particle_idx = threadIdx.x+blockIdx.x*blockDim.x;
    if(particle_idx >= pstruct.capacity)
        return;

    switch(pstruct.status_dev_p[particle_idx]) {
        case cuda_particle_struct::TERMINATED:
        case cuda_particle_struct::DETECTED:
        case cuda_particle_struct::PENDING:
            return;
        default:
            break;
    }

    const float3 pos = make_float3(
        pstruct.pos_x_dev_p[particle_idx],
        pstruct.pos_y_dev_p[particle_idx],
        pstruct.pos_z_dev_p[particle_idx]
    );

    if(!inside_AABB(pos, gstruct.AABB_center, gstruct.AABB_halfsize)) {
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
        return;
    }

    const int material_idx = pstruct.material_idx_dev_p[particle_idx];
    if(material_idx < 0) {
        // vacuum propgation
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::NO_EVENT;
        pstruct.distance_dev_p[particle_idx] =
            2.0f*norm3df(gstruct.AABB_halfsize.x, gstruct.AABB_halfsize.y, gstruct.AABB_halfsize.z);
        return;
    }

    const float K = pstruct.K_energy_dev_p[particle_idx];

    // material propagation
    if(K < mstruct.barrier_dev_p[material_idx]) {
        // EXIT: if energy is below barrier
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
        return;
    }

    float elastic_imfp, inelastic_imfp, total_imfp;
    {// log-log interpolate the elastic and inelastic inverse mean free path
        const float x = logf(K/mstruct.K_min)/logf(mstruct.K_max/mstruct.K_min);
        const int3 offset = make_int3(0, 0, material_idx);
        const int dim = mstruct.K_cnt;
        elastic_imfp = expf(interp1(mstruct.elastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, offset, dim, x));
        inelastic_imfp = expf(interp1(mstruct.inelastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, offset, dim, x));
        total_imfp = elastic_imfp+inelastic_imfp;
    }

    curand_wrapper curand(rand_state_dev_p[particle_idx]);

    pstruct.distance_dev_p[particle_idx] = -logf(curand.uniform())/total_imfp;
    if(curand.uniform() < elastic_imfp/total_imfp)
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::ELASTIC_EVENT;
    else
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::INELASTIC_EVENT;
}

__global__ void cuda_update_trajectory(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, scatter_options opt) {
    const int particle_idx = threadIdx.x+blockIdx.x*blockDim.x;
    if(particle_idx >= pstruct.capacity)
        return;

    switch(pstruct.status_dev_p[particle_idx]) {
        case cuda_particle_struct::TERMINATED:
        case cuda_particle_struct::DETECTED:
        case cuda_particle_struct::PENDING:
            return;
        default:
            break;
    }

    float3 pos = make_float3(
        pstruct.pos_x_dev_p[particle_idx],
        pstruct.pos_y_dev_p[particle_idx],
        pstruct.pos_z_dev_p[particle_idx]
    );

    float3 dir = make_float3(
        pstruct.dir_x_dev_p[particle_idx],
        pstruct.dir_y_dev_p[particle_idx],
        pstruct.dir_z_dev_p[particle_idx]
    );
    dir *= rnorm3df(dir.x, dir.y, dir.z);

    const int material_idx = pstruct.material_idx_dev_p[particle_idx];
    const int last_triangle_idx = pstruct.triangle_idx_dev_p[particle_idx];
    const float K = pstruct.K_energy_dev_p[particle_idx];

    float distance = pstruct.distance_dev_p[particle_idx];
    uint64_t location = 1;
    do {
        // define the root axis aligned bounding box
        float4 AABB = make_float4(
            gstruct.AABB_center.x,
            gstruct.AABB_center.y,
            gstruct.AABB_center.z,
            1.0f // scale factor
        );
        int index = 0;

        // traverse to location
        for(int i = 60-__clzll(location); i >= 0; i -= 3) {
            const unsigned int octant = (location>>i)&7;
            index = gstruct.octree_dev_p[index+octant];
            AABB.x += AABB.w*gstruct.AABB_halfsize.x*(1.0f*(octant&1)-0.5f);
            AABB.y += AABB.w*gstruct.AABB_halfsize.y*(1.0f*((octant&2)>>1)-0.5f);
            AABB.z += AABB.w*gstruct.AABB_halfsize.z*(1.0f*((octant&4)>>2)-0.5f);
            AABB.w *= 0.5f;
        }

        // traverse to leaf
        while(index >= 0) {
            unsigned int octant = 0;
            octant ^= (pos.x > AABB.x) ? 1 : 0;
            octant ^= (pos.y > AABB.y) ? 2 : 0;
            octant ^= (pos.z > AABB.z) ? 4 : 0;
            location = (location<<3)|octant;
            index = gstruct.octree_dev_p[index+octant];
            AABB.x += AABB.w*gstruct.AABB_halfsize.x*(1.0f*(octant&1)-0.5f);
            AABB.y += AABB.w*gstruct.AABB_halfsize.y*(1.0f*((octant&2)>>1)-0.5f);
            AABB.z += AABB.w*gstruct.AABB_halfsize.z*(1.0f*((octant&4)>>2)-0.5f);
            AABB.w *= 0.5f;
        }
        index = -index;

        // determine cell/triangle intersections for current leaf
        float intersect;
        int target; // target >= 0 : triangle index
                    // target = -1 : cell x intersection
                    // target = -2 : cell y intersection
                    // target = -4 : cell z intersection

        const float3 t = AABB_intersect(pos, dir, make_float3(AABB.x, AABB.y, AABB.z), AABB.w*gstruct.AABB_halfsize);
        if((t.x < t.y) && (t.x < t.z)) {
            intersect = t.x;
            target = -1;
        } else if((t.y < t.x) && (t.y < t.z)) {
            intersect = t.y;
            target = -2;
        } else {
            intersect = t.z;
            target = -4;
        }

        while(true) {
            const int triangle_idx = gstruct.octree_dev_p[index++];
            if(triangle_idx < 0)
                break;
            if(triangle_idx == last_triangle_idx)
                continue;

            const float3 O = make_float3(
                gstruct.triangle_r0x_dev_p[triangle_idx],
                gstruct.triangle_r0y_dev_p[triangle_idx],
                gstruct.triangle_r0z_dev_p[triangle_idx]
            );
            const float3 e1 = make_float3(
                gstruct.triangle_e1x_dev_p[triangle_idx],
                gstruct.triangle_e1y_dev_p[triangle_idx],
                gstruct.triangle_e1z_dev_p[triangle_idx]
            );
            const float3 e2 = make_float3(
                gstruct.triangle_e2x_dev_p[triangle_idx],
                gstruct.triangle_e2y_dev_p[triangle_idx],
                gstruct.triangle_e2z_dev_p[triangle_idx]
            );

            int mat_idx_out;
            if(dot_product(cross_product(e1, e2), dir) < 0)
                mat_idx_out = gstruct.material_idx_in_dev_p[triangle_idx];
            else
                mat_idx_out = gstruct.material_idx_out_dev_p[triangle_idx];

            if((mat_idx_out == material_idx) || (mat_idx_out == triangle::NOP))
                continue;
            else if((mat_idx_out == triangle::DETECTOR_LT50) && (K >= 50))
                continue;
            else if((mat_idx_out == triangle::DETECTOR_GE50) && (K < 50))
                continue;

            const float t = triangle_intersect(pos, dir, O, e1, e2);
            if((t > 0) && (t <= intersect+_eps)) {
                intersect = t;
                target = triangle_idx;
            }
        }

        // manage calculated intersections
        if(intersect >= distance) {
            // EXIT: distance traveled without triangle intersections
            pstruct.pos_x_dev_p[particle_idx] = pos.x+dir.x*distance;
            pstruct.pos_y_dev_p[particle_idx] = pos.y+dir.y*distance;
            pstruct.pos_z_dev_p[particle_idx] = pos.z+dir.z*distance;
            return;
        } else if(target >= 0) {
            // EXIT: triangle intersection
            pstruct.status_dev_p[particle_idx] = cuda_particle_struct::INTERSECT_EVENT;
            pstruct.triangle_idx_dev_p[particle_idx] = target;
            pstruct.pos_x_dev_p[particle_idx] = pos.x+dir.x*intersect;
            pstruct.pos_y_dev_p[particle_idx] = pos.y+dir.y*intersect;
            pstruct.pos_z_dev_p[particle_idx] = pos.z+dir.z*intersect;
            return;
        }
        distance -= intersect;
        pos.x += dir.x*intersect;
        pos.y += dir.y*intersect;
        pos.z += dir.z*intersect;

        // find adjacent node
        unsigned int mask = -target;
        unsigned int value;
        // mask = 1 (001), value = 0 (000) : xx0 --> xx1 (find neighbor in positive x direction)
        // mask = 1 (001), value = 1 (000) : xx1 --> xx0 (find neighbor in negative x direction)
        // mask = 2 (010), value = 0 (000) : x0x --> x1x (find neighbor in positive y direction)
        // mask = 2 (010), value = 2 (010) : x1x --> x0x (find neighbor in negative y direction)
        // mask = 4 (100), value = 0 (000) : 0xx --> 1xx (find neighbor in positive z direction)
        // mask = 4 (100), value = 4 (100) : 1xx --> 0xx (find neighbor in negative z direction)
        if(mask == 1)
            value = (dir.x >= 0) ? 0 : 1;
        else if(mask == 2)
            value = (dir.y >= 0) ? 0 : 2;
        else
            value = (dir.z >= 0) ? 0 : 4;
        while(location > 1) {
            if((location&mask) == value) {
                location ^= mask;
                break;
            }
            location >>= 3;
        }

    } while(location > 1);

    // EXIT: particle is out of grid
    pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
    return;
}

__global__ void cuda_intersection_event(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p, scatter_options opt) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= pstruct.capacity)
        return;
    const int particle_idx = pstruct.particle_idx_dev_p[i];

    switch(pstruct.status_dev_p[particle_idx]) {
        case cuda_particle_struct::INTERSECT_EVENT:
            break;
        default:
            return;
    }

    float3 dir = make_float3(
        pstruct.dir_x_dev_p[particle_idx],
        pstruct.dir_y_dev_p[particle_idx],
        pstruct.dir_z_dev_p[particle_idx]
    );
    dir *= rnorm3df(dir.x, dir.y, dir.z);

    const int triangle_idx = pstruct.triangle_idx_dev_p[particle_idx];
    const float3 e1 = make_float3(
        gstruct.triangle_e1x_dev_p[triangle_idx],
        gstruct.triangle_e1y_dev_p[triangle_idx],
        gstruct.triangle_e1z_dev_p[triangle_idx]
    );
    const float3 e2 = make_float3(
        gstruct.triangle_e2x_dev_p[triangle_idx],
        gstruct.triangle_e2y_dev_p[triangle_idx],
        gstruct.triangle_e2z_dev_p[triangle_idx]
    );
    float3 normal = cross_product(e1, e2);
    normal *= rnorm3df(normal.x, normal.y, normal.z);
    const float cos_theta = dot_product(normal, dir);
    const float sin_theta = sqrtf(1.0f-__saturatef(cos_theta*cos_theta));

    int material_idx_in, material_idx_out;
    if(cos_theta > 0) {
        material_idx_in = gstruct.material_idx_in_dev_p[triangle_idx];
        material_idx_out = gstruct.material_idx_out_dev_p[triangle_idx];
    } else {
        material_idx_in = gstruct.material_idx_out_dev_p[triangle_idx];
        material_idx_out = gstruct.material_idx_in_dev_p[triangle_idx];
    }

    switch(material_idx_out) {
        case triangle::DETECTOR:
        case triangle::DETECTOR_LT50:
        case triangle::DETECTOR_GE50:
            pstruct.status_dev_p[particle_idx] = cuda_particle_struct::DETECTED;
            return;
        case triangle::TERMINATOR:
            pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
            return;
        case triangle::MIRROR:
            pstruct.dir_x_dev_p[particle_idx] = dir.x-2.0f*normal.x*cos_theta;
            pstruct.dir_y_dev_p[particle_idx] = dir.y-2.0f*normal.y*cos_theta;
            pstruct.dir_z_dev_p[particle_idx] = dir.z-2.0f*normal.z*cos_theta;
            return;
        default:
            break;
    }

    float dU = 0;
    float band_edge[2] = {0.0f, 0.0f};
    float effective_mass[2] = {1.0f, 1.0f};
    if(material_idx_in >= 0) {
        dU -= mstruct.barrier_dev_p[material_idx_in];
        band_edge[0] = mstruct.band_edge_dev_p[material_idx_in];
        effective_mass[0] = mstruct.effective_mass_dev_p[material_idx_in];
    }
    if(material_idx_out >= 0) {
        dU += mstruct.barrier_dev_p[material_idx_out];
        band_edge[1] = mstruct.band_edge_dev_p[material_idx_out];
        effective_mass[1] = mstruct.effective_mass_dev_p[material_idx_out];
    }

    curand_wrapper curand(rand_state_dev_p[particle_idx]);

    const float K = pstruct.K_energy_dev_p[particle_idx];
    if(K*cos_theta*cos_theta+dU > 0) {
        const float r = effective_mass[1]/effective_mass[0];
        const float s = sqrtf(r*(K+dU-band_edge[1])/(K-band_edge[0])-sin_theta*sin_theta)/fabsf(cos_theta);
        const float T = (opt.quantum_transmission_flag) ? 4.0f*s/((1.0f+s)*(1.0f+s)) : 1.0f;
        if(curand.uniform() < T) {
            if(opt.interface_refraction_flag) {
                pstruct.dir_x_dev_p[particle_idx] = (dir.x-normal.x*cos_theta)+s*normal.x*cos_theta;
                pstruct.dir_y_dev_p[particle_idx] = (dir.y-normal.y*cos_theta)+s*normal.y*cos_theta;
                pstruct.dir_z_dev_p[particle_idx] = (dir.z-normal.z*cos_theta)+s*normal.z*cos_theta;
            }
            pstruct.K_energy_dev_p[particle_idx] = K+dU;
            pstruct.material_idx_dev_p[particle_idx] = material_idx_out;
            return;
        }
    }
    if((opt.interface_absorption_flag) && (dU < 0) && (curand.uniform() < expf(1.0f+0.5f*K/dU))) {
        // surface absorption? (this is in accordance with Kieft & Bosch code)
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
        return;
    }
    // total internal reflection
    pstruct.dir_x_dev_p[particle_idx] = dir.x-2.0f*normal.x*cos_theta;
    pstruct.dir_y_dev_p[particle_idx] = dir.y-2.0f*normal.y*cos_theta;
    pstruct.dir_z_dev_p[particle_idx] = dir.z-2.0f*normal.z*cos_theta;
}

__global__ void cuda_elastic_event(cuda_particle_struct pstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p, scatter_options opt) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= pstruct.capacity)
        return;
    const int particle_idx = pstruct.particle_idx_dev_p[i];

    switch(pstruct.status_dev_p[particle_idx]) {
        case cuda_particle_struct::ELASTIC_EVENT:
            break;
        default:
            return;
    }

    // forget last intersected triangle
    pstruct.triangle_idx_dev_p[particle_idx] = -1;

    curand_wrapper curand(rand_state_dev_p[particle_idx]);

    const int material_idx = pstruct.material_idx_dev_p[particle_idx];
    const float K = pstruct.K_energy_dev_p[particle_idx];

    float cos_theta, sin_theta;
    {// sample random elastic scatter angle
        const float x = logf(K/mstruct.K_min)/logf(mstruct.K_max/mstruct.K_min);
        const float y = curand.uniform();
        const int3 offset = make_int3(0, 1, material_idx);
        const int2 dim = make_int2(mstruct.K_cnt, mstruct.P_cnt);
        cos_theta = interp2(mstruct.elastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, offset, dim, x, y);
        cos_theta = clamp(cos_theta, -1.0f, 1.0f);
        sin_theta = sqrtf(1.0f-cos_theta*cos_theta);
    }

    float3 dir = make_float3(
        pstruct.dir_x_dev_p[particle_idx],
        pstruct.dir_y_dev_p[particle_idx],
        pstruct.dir_z_dev_p[particle_idx]
    );
    dir *= rnorm3df(dir.x, dir.y, dir.z);

    float3 normal_dir = make_normal_vec(dir, 2.0f*_pi*curand.uniform());
    normal_dir *= rnorm3df(normal_dir.x, normal_dir.y, normal_dir.z);

    pstruct.dir_x_dev_p[particle_idx] = dir.x*cos_theta+normal_dir.x*sin_theta;
    pstruct.dir_y_dev_p[particle_idx] = dir.y*cos_theta+normal_dir.y*sin_theta;
    pstruct.dir_z_dev_p[particle_idx] = dir.z*cos_theta+normal_dir.z*sin_theta;

    if(K < 200) {
        if(opt.acoustic_phonon_loss_flag)
            pstruct.K_energy_dev_p[particle_idx] = K-fminf(50e-3f, mstruct.phonon_loss_dev_p[material_idx]);
    } else {
        if(opt.atomic_recoil_loss_flag) {
            #warning "fixed energy loss due to atom recoil assumes silicon material"
            pstruct.K_energy_dev_p[particle_idx] = K-2.0f*(_me*_NA)*K*(1.0f-cos_theta)/28.1f;
        }
    }
}

__global__ void cuda_inelastic_event(cuda_particle_struct pstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p, scatter_options opt) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= pstruct.capacity)
        return;
    const int primary_idx = pstruct.particle_idx_dev_p[i];

    switch(pstruct.status_dev_p[primary_idx]) {
        case cuda_particle_struct::INELASTIC_EVENT:
        case cuda_particle_struct::PENDING:
            break;
        default:
            return;
    }

    const int secondary_idx = pstruct.particle_idx_dev_p[pstruct.capacity-1-i];
    if(pstruct.status_dev_p[secondary_idx] != cuda_particle_struct::TERMINATED) {
        pstruct.status_dev_p[primary_idx] = cuda_particle_struct::PENDING;
        return;
    }
    pstruct.status_dev_p[primary_idx] = cuda_particle_struct::INELASTIC_EVENT;

    // forget last intersected triangle
    pstruct.triangle_idx_dev_p[primary_idx] = -1;

    const int material_idx = pstruct.material_idx_dev_p[primary_idx];
    const float fermi = mstruct.fermi_dev_p[material_idx];
    const float K = pstruct.K_energy_dev_p[primary_idx];

    curand_wrapper curand(rand_state_dev_p[primary_idx]);

    float omega0;
    {// sample random zero-momentum energy loss of the primary electron
        const float x = logf(K/mstruct.K_min)/logf(mstruct.K_max/mstruct.K_min);
        const float y = curand.uniform();
        const int3 offset = make_int3(0, 1, material_idx);
        const int2 dim = make_int2(mstruct.K_cnt, mstruct.P_cnt);
        omega0 = expf(interp2(mstruct.inelastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, offset, dim, x, y));
        omega0 = posf(omega0);
    }

    float binding;
    {// sample random binding energy of the secondary electron
        const float x = logf(omega0/mstruct.K_min)/logf(mstruct.K_max/mstruct.K_min);
        const float y = curand.uniform();
        const int ix = clamp(__float2int_rd(x*(mstruct.K_cnt-1)), 0, mstruct.K_cnt-1);
        const int iy = clamp(__float2int_rd(y*(mstruct.P_cnt-1)), 0, mstruct.P_cnt-1);
        binding = cuda_make_ptr(mstruct.ionization_dev_p, mstruct.pitch, mstruct.P_cnt+1, iy, material_idx)[ix];
    }

    float omega;
    {// sample random total energy loss for the primary electron
        float omega_max = 0.5f*(K+omega0-fermi); // upper limit of eq. 9 in Ashley, but corrected for the fermi energy
        float omega_min = omega0;
        float w0 = fminf(omega0-1.0f, posf(binding)-fermi);
        if(K > 2.0f*omega0) {
            omega_min = 0.5f*(K+omega0-sqrtf(K*(K-2.0f*omega0))); // equation 10 in Ashley
            w0 = omega0;
        }
        const float U = curand.uniform();
        if((w0 > 0) && (omega_min > w0) && (omega_min < omega_max)) {
            // For nonzero binding energy, sample omega according to equation 7 in Ashley,
            // using the lower and upper limits as defined above.
            // For inner-shell ionization (Ebind > 50 eV) we substitute the Fermi-energy corrected
            // binding energy for omegaprime (so that the differential cross section becomes inversely
            // proportional to both the total energy transfer and the kinetic energy of the secondary
            // electron).
            const float f_min = 1.0f/w0*logf((omega_min-w0)/omega_min);
            const float f_max = 1.0f/w0*logf((omega_max-w0)/omega_max);
            omega = -w0/expm1f(w0*(f_min*(1.0f-U)+f_max*U));
        } else {
            // In some cases (typically only occuring for binding < 50 eV) we get omega_min > omega_max.
            // This is due to our Fermi energy correction in the definition of omega_max. Physically, this
            // means that momentum cannot be conserved because the primary electron cannot have a final
            // kinetic energy that is lower than the Fermi energy. In this (relatively rare) case we have
            // to ignore momentum conservation and probe omega according to a 1/(omega)^2 distribution
            // with omega0 and omega_max as lower and upper limits respectively.
            omega = omega0*omega_max/(omega0*(1.0f-U)+omega_max*U);
        }
    }

    if(binding < 0) {
        const float band_gap = mstruct.band_gap_dev_p[material_idx];
        if(band_gap < 0) {
            // metal: excitation of a fermi sea electron
            // TODO
        } else if(omega0 > band_gap) {
            // electron excitation across the band gap
            binding = band_gap;
        } else {
            // sub-band gap energy loss in semiconductors and insulators
            // energy loss due to longitudinal optical phonon excitation is assumed
            if(opt.optical_phonon_loss_flag)
                pstruct.K_energy_dev_p[primary_idx] = K-omega0;
            return;
        }
    }
    binding = posf(binding);

    float3 primary_dir = make_float3(
        pstruct.dir_x_dev_p[primary_idx],
        pstruct.dir_y_dev_p[primary_idx],
        pstruct.dir_z_dev_p[primary_idx]
    );
    primary_dir *= rnorm3df(primary_dir.x, primary_dir.y, primary_dir.z);

    float3 normal_dir = make_normal_vec(primary_dir, 2.0f*_pi*curand.uniform());
    normal_dir *= rnorm3df(normal_dir.x, normal_dir.y, normal_dir.z);

    const float _K = K-fermi+2.0f*binding;
    const float dK = binding+omega;
    const float cos_theta = sqrtf(dK/_K);
    const float sin_theta = sqrtf(1.0f-__saturatef(cos_theta*cos_theta));

    float3 secondary_dir = make_float3(
        primary_dir.x*cos_theta+normal_dir.x*sin_theta,
        primary_dir.y*cos_theta+normal_dir.y*sin_theta,
        primary_dir.z*cos_theta+normal_dir.z*sin_theta
    );
    if(opt.instantaneous_momentum_flag) {
        const float random_cos = 2.0f*curand.uniform()-1.0f;
        const float random_phi = 2.0f*_pi*curand.uniform();
        secondary_dir *= rnorm3df(secondary_dir.x, secondary_dir.y, secondary_dir.z);
        secondary_dir += sqrtf(binding/dK)*make_unit_vec(random_cos, random_phi);
    }
    secondary_dir *= rnorm3df(secondary_dir.x, secondary_dir.y, secondary_dir.z);

    if(opt.generate_secondary_flag) {
        pstruct.status_dev_p[secondary_idx] = cuda_particle_struct::NEW_SECONDARY;
        pstruct.material_idx_dev_p[secondary_idx] = material_idx;
        pstruct.particle_tag_dev_p[secondary_idx] = pstruct.particle_tag_dev_p[primary_idx];

        pstruct.pos_x_dev_p[secondary_idx] = pstruct.pos_x_dev_p[primary_idx];
        pstruct.pos_y_dev_p[secondary_idx] = pstruct.pos_y_dev_p[primary_idx];
        pstruct.pos_z_dev_p[secondary_idx] = pstruct.pos_z_dev_p[primary_idx];

        pstruct.K_energy_dev_p[secondary_idx] = fermi+omega-binding;
        pstruct.dir_x_dev_p[secondary_idx] = secondary_dir.x;
        pstruct.dir_y_dev_p[secondary_idx] = secondary_dir.y;
        pstruct.dir_z_dev_p[secondary_idx] = secondary_dir.z;
    }

    // primary direction determined by non-relativistic momentum-conservation, i.e.:
    //   sin(theta)*primary_dir_2 = primary_dir - cos(theta)*secondary_dir;
    pstruct.K_energy_dev_p[primary_idx] = K-omega;
    if(opt.momentum_conservation_flag) {
        pstruct.dir_x_dev_p[primary_idx] = primary_dir.x-cos_theta*secondary_dir.x;
        pstruct.dir_y_dev_p[primary_idx] = primary_dir.y-cos_theta*secondary_dir.y;
        pstruct.dir_z_dev_p[primary_idx] = primary_dir.z-cos_theta*secondary_dir.z;
    }
}