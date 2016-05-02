/**
 * @file src/cdsem/cuda_kernels.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "cuda_kernels.cuh"
#include <cfloat>
#include <common/cuda_make_ptr.cuh>
#include <common/cuda_vec3_math.cuh>

__device__ const float eps = 10.0f*FLT_EPSILON;
__device__ const float mc2 = 5.1099897e+05f;
__device__ const float pi = 6.2831853e+00f;

/*
namespace fermi_sea {
__device__ const int n = 32;
__device__ const float z1 = 1.0000000e-03f;
__device__ const float z2 = 1.0000000e+01f;
__device__ const float p1 = 2.1088175e-05f;
__device__ const float p2 = 5.4595325e+01f;
__device__ const float z2p[n] = {
    -1.0766798e+01f, -1.0321033e+01f, -9.8752317e+00f, -9.4293823e+00f,
    -8.9834681e+00f, -8.5374670e+00f, -8.0913477e+00f, -7.6450720e+00f,
    -7.1985850e+00f, -6.7518148e+00f, -6.3046651e+00f, -5.8570080e+00f,
    -5.4086747e+00f, -4.9594412e+00f, -4.5090156e+00f, -4.0570192e+00f,
    -3.6029665e+00f, -3.1462471e+00f, -2.6861079e+00f, -2.2216477e+00f,
    -1.7518265e+00f, -1.2755001e+00f, -7.9148793e-01f, -2.9867163e-01f,
    +2.0388594e-01f, +7.1682256e-01f, +1.2403889e+00f, +1.7744150e+00f,
    +2.3183398e+00f, +2.8712943e+00f, +3.4322128e+00f, +3.9999483e+00f
};
__device__ const float p2z[n] = {
    -6.9077554e+00f, -6.5902653e+00f, -6.2728038e+00f, -5.9553800e+00f,
    -5.6380086e+00f, -5.3207092e+00f, -5.0035086e+00f, -4.6864424e+00f,
    -4.3695607e+00f, -4.0529299e+00f, -3.7366409e+00f, -3.4208145e+00f,
    -3.1056123e+00f, -2.7912462e+00f, -2.4779902e+00f, -2.1661916e+00f,
    -1.8562782e+00f, -1.5487578e+00f, -1.2442060e+00f, -9.4323510e-01f,
    -6.4644718e-01f, -3.5437316e-01f, -6.7411400e-02f, +2.1421888e-01f,
    +4.9049795e-01f, +7.6159561e-01f, +1.0278367e+00f, +1.2896539e+00f,
    +1.5475388e+00f, +1.8019993e+00f, +2.0535290e+00f, +2.3025851e+00f
};
};
*/

template<typename T>
__device__ T clamp(T x, T x1, T x2) {
    return max(x1, min(x2, x));
}

__device__ float interp1(const float* ptr, int pitch, int height, int ix, int iy, int iz, float sx) {
    ptr = cuda_make_ptr(ptr, pitch, height, iy, iz);
    return (1.0f-sx)*ptr[ix]+sx*ptr[ix+1];
}

__device__ float interp2(const float* ptr, int pitch, int height, int ix, int iy, int iz, float sx, float sy) {
    return sy*interp1(ptr, pitch, height, ix, iy, iz, sx)+(1.0f-sy)*interp1(ptr, pitch, height, ix, iy+1, iz, sx);
}

__device__ bool inside_AABB(float3 pos, float3 center, float3 halfsize) {
    if((pos.x > center.x-halfsize.x) && (pos.x < center.x+halfsize.x))
    if((pos.y > center.y-halfsize.y) && (pos.y < center.y+halfsize.y))
    if((pos.z > center.z-halfsize.z) && (pos.z < center.z+halfsize.z))
        return true;
    return false;
}

__global__ void cuda_init_rand_state(curandState* rand_state_p, unsigned long long seed, int n) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= n)
        return;
    curand_init(seed, i, 0, &(rand_state_p[i]));
}

__global__ void cuda_init_trajectory(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p) {
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

    // material propagation
    if(pstruct.K_energy_dev_p[particle_idx] < mstruct.barrier_dev_p[material_idx]) {
        // EXIT: if energy is below barrier
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
        return;
    }

    curandState rand_state = rand_state_dev_p[particle_idx];

    const float K = pstruct.K_energy_dev_p[particle_idx];
    const float x = __fdividef(__logf(K/mstruct.K_min), __logf(mstruct.K_max/mstruct.K_min))*(mstruct.K_cnt-1);
    const int ix = clamp(__float2int_rd(x), 0, mstruct.K_cnt-2);

    const float elastic_imfp =
        __expf(interp1(mstruct.elastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, ix, 0, material_idx, x-ix));
    const float inelastic_imfp =
        __expf(interp1(mstruct.inelastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, ix, 0, material_idx, x-ix));
    const float total_imfp = elastic_imfp+inelastic_imfp;
    pstruct.distance_dev_p[particle_idx] = __fdividef(-__logf(curand_uniform(&rand_state)), total_imfp);

    if(curand_uniform(&rand_state) < elastic_imfp/total_imfp)
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::ELASTIC_EVENT;
    else
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::INELASTIC_EVENT;

    rand_state_dev_p[particle_idx] = rand_state;
}

__global__ void cuda_update_trajectory(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct) {
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
    dir = dir*rnorm3df(dir.x, dir.y, dir.z);

    const int material_idx = pstruct.material_idx_dev_p[particle_idx];
    const int last_triangle_idx = pstruct.triangle_idx_dev_p[particle_idx];
    const float K = pstruct.K_energy_dev_p[particle_idx];

    float distance = pstruct.distance_dev_p[particle_idx];
    uint64_t location = 1;
    do {
        // traverse to location
        float4 AABB = make_float4(
            gstruct.AABB_center.x,
            gstruct.AABB_center.y,
            gstruct.AABB_center.z,
            1.0f // scale factor
        );
        int index = 0;
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
                    // target = -1 : node x intersection
                    // target = -2 : node y intersection
                    // target = -4 : node z intersection

        const float3 t = make_float3(
            __fdividef(AABB.x+copysignf(AABB.w*gstruct.AABB_halfsize.x+eps, dir.x)-pos.x, dir.x),
            __fdividef(AABB.y+copysignf(AABB.w*gstruct.AABB_halfsize.y+eps, dir.y)-pos.y, dir.y),
            __fdividef(AABB.z+copysignf(AABB.w*gstruct.AABB_halfsize.z+eps, dir.z)-pos.z, dir.z)
        );
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

            // T. Möller and B. Trumbore, Journal of Graphics Tools, 2(1):21--28, 1997.
            const float3 pvec = cross_product(dir, e2);
            const float det = dot_product(e1, pvec);
            if(fabsf(det) < eps)
                continue;
            const float3 tvec = pos-make_float3(
                gstruct.triangle_r0x_dev_p[triangle_idx],
                gstruct.triangle_r0y_dev_p[triangle_idx],
                gstruct.triangle_r0z_dev_p[triangle_idx]
            );
            const float u = __fdividef(dot_product(tvec, pvec), det);
            if((u < -eps) || (u > 1.0f+eps))
                continue;
            const float3 qvec = cross_product(tvec, e1);
            const float v = __fdividef(dot_product(dir, qvec), det);
            if((v < -eps) || (u+v > 1.0f+eps))
                continue;
            const float t = __fdividef(dot_product(e2, qvec), det);
            if((t > 0) && (t <= intersect+eps)) {
                intersect = t;
                target = triangle_idx;
            }
        }

        // manage calculated intersections
        if(intersect >= distance) {
            // EXIT: distance travelled without triangle intersections
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

__global__ void cuda_intersection_event(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p) {
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
    dir = dir*rnorm3df(dir.x, dir.y, dir.z);

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
    normal = normal*rnorm3df(normal.x, normal.y, normal.z);
    const float cos_alpha = dot_product(normal, dir);

    int material_idx_in, material_idx_out;
    if(cos_alpha > 0) {
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
            pstruct.dir_x_dev_p[particle_idx] = dir.x-2.0f*normal.x*cos_alpha;
            pstruct.dir_y_dev_p[particle_idx] = dir.y-2.0f*normal.y*cos_alpha;
            pstruct.dir_z_dev_p[particle_idx] = dir.z-2.0f*normal.z*cos_alpha;
            return;
        default:
            break;
    }

    float dU = 0;
    if(material_idx_out >= 0)
        dU += mstruct.barrier_dev_p[material_idx_out];
    if(material_idx_in >= 0)
        dU -= mstruct.barrier_dev_p[material_idx_in];

    curandState rand_state = rand_state_dev_p[particle_idx];

    // R. Shimizu and Z. J. Ding, Rep. Prog. Phys., 55, 487-531, 1992
    //  see equations 3.20, 3.23 and 3.24
    const float K = pstruct.K_energy_dev_p[particle_idx];
    const float z = sqrtf(fmaxf(0.0f, 1.0f+dU/(K*cos_alpha*cos_alpha)));
    if((K*cos_alpha*cos_alpha+dU > 0) && (curand_uniform(&rand_state) < __fdividef(4.0f*z, ((1.0f+z)*(1.0f+z))))) {
        pstruct.dir_x_dev_p[particle_idx] = (dir.x-normal.x*cos_alpha)+normal.x*cos_alpha*z;
        pstruct.dir_y_dev_p[particle_idx] = (dir.y-normal.y*cos_alpha)+normal.y*cos_alpha*z;
        pstruct.dir_z_dev_p[particle_idx] = (dir.z-normal.z*cos_alpha)+normal.z*cos_alpha*z;
        pstruct.K_energy_dev_p[particle_idx] = K+dU;
        pstruct.material_idx_dev_p[particle_idx] = material_idx_out;
    } else if((dU < 0) && (curand_uniform(&rand_state) < __expf(1.0f+0.5f*K/dU))) {
        // surface absorption? (this is in accordance with Kieft & Bosch code)
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
    } else {
        // total internal reflection
        pstruct.dir_x_dev_p[particle_idx] = dir.x-2.0f*normal.x*cos_alpha;
        pstruct.dir_y_dev_p[particle_idx] = dir.y-2.0f*normal.y*cos_alpha;
        pstruct.dir_z_dev_p[particle_idx] = dir.z-2.0f*normal.z*cos_alpha;
    }

    rand_state_dev_p[particle_idx] = rand_state;
}

__global__ void cuda_elastic_event(cuda_particle_struct pstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p) {
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

    curandState rand_state = rand_state_dev_p[particle_idx];

    const int material_idx = pstruct.material_idx_dev_p[particle_idx];
    const float K = pstruct.K_energy_dev_p[particle_idx];
    const float x = __fdividef(__logf(K/mstruct.K_min), __logf(mstruct.K_max/mstruct.K_min))*(mstruct.K_cnt-1);
    const int ix = clamp(__float2int_rd(x), 0, mstruct.K_cnt-2);
    const float y = curand_uniform(&rand_state)*(mstruct.P_cnt-1);
    const int iy = clamp(__float2int_rd(y), 0, mstruct.P_cnt-2);

    const float cos_theta = clamp(interp2(mstruct.elastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, ix, iy+1, material_idx, x-ix, y-iy), -1.0f, 1.0f);
    const float sin_theta = sqrtf(1.0f-cos_theta*cos_theta);

    float3 dir = make_float3(
        pstruct.dir_x_dev_p[particle_idx],
        pstruct.dir_y_dev_p[particle_idx],
        pstruct.dir_z_dev_p[particle_idx]
    );
    dir = dir*rnorm3df(dir.x, dir.y, dir.z);

    float sin_azimuth, cos_azimuth;
    __sincosf(atan2f(dir.y, dir.x), &sin_azimuth, &cos_azimuth);

    float sin_phi, cos_phi;
    __sincosf(2.0f*pi*curand_uniform(&rand_state), &sin_phi, &cos_phi);

    const float3 unit_v = make_float3(
        dir.z*cos_azimuth,
        dir.z*sin_azimuth,
        -sqrtf(fmaxf(0.0f, 1.0f-dir.z*dir.z))
    );
    const float3 unit_u = cross_product(unit_v, dir);

    pstruct.dir_x_dev_p[particle_idx] = dir.x*cos_theta+(unit_u.x*cos_phi+unit_v.x*sin_phi)*sin_theta;
    pstruct.dir_y_dev_p[particle_idx] = dir.y*cos_theta+(unit_u.y*cos_phi+unit_v.y*sin_phi)*sin_theta;
    pstruct.dir_z_dev_p[particle_idx] = dir.z*cos_theta+(unit_u.z*cos_phi+unit_v.z*sin_phi)*sin_theta;

    rand_state_dev_p[particle_idx] = rand_state;
}

__global__ void cuda_inelastic_event(cuda_particle_struct pstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p) {
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

    curandState rand_state = rand_state_dev_p[primary_idx];

    const int material_idx = pstruct.material_idx_dev_p[primary_idx];
    const float fermi = mstruct.fermi_dev_p[material_idx];
    const float K = pstruct.K_energy_dev_p[primary_idx];
    float x, y;
    int ix, iy;

    x = __fdividef(__logf(K/mstruct.K_min), __logf(mstruct.K_max/mstruct.K_min))*(mstruct.K_cnt-1);
    ix = clamp(__float2int_rd(x), 0, mstruct.K_cnt-2);
    y = curand_uniform(&rand_state)*(mstruct.P_cnt-1);
    iy = clamp(__float2int_rd(y), 0, mstruct.P_cnt-2);
    const float omega0 = clamp(interp2(mstruct.inelastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, ix, iy+1, material_idx, x-ix, y-iy), 0.0f, K-fermi);

    x = __fdividef(__logf(omega0/mstruct.K_min), __logf(mstruct.K_max/mstruct.K_min))*(mstruct.K_cnt-1);
    ix = clamp(__float2int_rd(x), 0, mstruct.K_cnt-1);
    y = curand_uniform(&rand_state)*(mstruct.P_cnt-1);
    iy = __float2int_rd(y);
    float binding = cuda_make_ptr(mstruct.ionization_dev_p, mstruct.pitch, mstruct.P_cnt+1, iy+1, material_idx)[ix];

    float omega;
    // upper limit of equation 9 in Ashley, but corrected for the fermi energy.
    float omega_max = 0.5f*(K+omega0-fermi);
    float omega_min = omega0;
    float w0 = min(omega0-1.0f, fmaxf(0.0f, binding)-fermi);
    if(K > 2.0f*omega0) {
        // equation 10 in Ashley
        omega_min = 0.5f*K*(1.0f-sqrtf(1.0f-2.0f*omega0/K)+omega0/K);
        w0 = omega0;
    }
    if((w0 > 0) && (omega_min > w0) && (omega_min < omega_max)) {
        // For nonzero binding energy, sample omega according to equation 7 in Ashley,
        // using the lower and upper limits as defined above.
        // For inner-shell ionization (Ebind > 50 eV) we substitute the Fermi-energy corrected
        // binding energy for omegaprime (so that the differential cross section becomes inversely
        // proportional to both the total energy transfer and the kinetic energy of the secondary
        // electron).
        const float U = curand_uniform(&rand_state);
        omega = w0/(1.0f-(1.0f-w0/omega_min)*__expf(U*log1pf(-w0/omega_max))*__expf(-U*log1pf(-w0/omega_min)));
    } else {
        // In some cases (typically only occuring for binding < 50 eV) we get omega_min > omega_max.
        // This is due to our Fermi energy correction in the definition of omega_max. Physically, this
        // means that momentum cannot be conserved because the primary electron cannot have a final
        // kinetic energy that is lower than the Fermi energy. In this (relatively rare) case we have
        // to ignore momentum conservation and probe omega according to a 1/(omega)^2 distribution
        // with omega0 and omega_max as lower and upper limits respectively.
        const float U = curand_uniform(&rand_state);
        omega = omega0/(1.0f-U*(1.0f-omega0/omega_max));
    }

    if(binding < 0) {
        const float bandgap = mstruct.bandgap_dev_p[material_idx];
        if(bandgap < 0) {
            // metal: excitation of a fermi sea electron
            /*
            const float z = F/omega;
            if((z-1.0f > fermi_sea::z1) && (z < fermi_sea::z2)) {
                const float x1 = logscale(z-1.0f,fermi_sea::z1,fermi_sea::z2);
                const float x2 = logscale(z,fermi_sea::z1,fermi_sea::z2);
                const float p1 = __expf(interp1(fermi_sea::z2p,fermi_sea::n,x1));
                const float p2 = __expf(interp1(fermi_sea::z2p,fermi_sea::n,x2));
                const float U = curand_uniform(&rand_state);
                const float y = logscale(p1*(1.0f-U)+p2*U,fermi_sea::p1,fermi_sea::p2);
                binding = F-omega*__expf(interp1(fermi_sea::p2z,fermi_sea::n,y));
            }
            */
        } else if(omega0 > bandgap) {
            // electron excitation across the bandgap
            binding = bandgap;
        } else {
            // sub-bandgap energy loss in semiconductors and insulators
            // energy loss due to longitudinal optical phonon excitation is assumed
            pstruct.K_energy_dev_p[primary_idx] = K-omega0;
            rand_state_dev_p[primary_idx] = rand_state;
            return;
        }
    }

    const float _K = K-fermi+2.0f*binding;
    const float dK = binding+omega;
    const float cos_alpha = sqrtf(__saturatef(__fdividef((1.0f-dK/_K)*(1.0f+0.5f*_K/mc2), 1.0f+0.5f*(_K-dK)/mc2)));
    const float sin_alpha = sqrtf(1.0f-cos_alpha*cos_alpha);
    const float cos_beta = sqrtf(__saturatef(__fdividef((dK/_K)*(1.0f+0.5f*_K/mc2), 1.0f+0.5f*dK/mc2)));
    const float sin_beta = sqrtf(1.0f-cos_beta*cos_beta);

    pstruct.status_dev_p[secondary_idx] = cuda_particle_struct::NEW_SECONDARY;
    pstruct.material_idx_dev_p[secondary_idx] = material_idx;
    pstruct.particle_tag_dev_p[secondary_idx] = pstruct.particle_tag_dev_p[primary_idx];
    pstruct.pos_x_dev_p[secondary_idx] = pstruct.pos_x_dev_p[primary_idx];
    pstruct.pos_y_dev_p[secondary_idx] = pstruct.pos_y_dev_p[primary_idx];
    pstruct.pos_z_dev_p[secondary_idx] = pstruct.pos_z_dev_p[primary_idx];

    float3 dir = make_float3(
        pstruct.dir_x_dev_p[primary_idx],
        pstruct.dir_y_dev_p[primary_idx],
        pstruct.dir_z_dev_p[primary_idx]
    );
    dir = dir*rnorm3df(dir.x, dir.y, dir.z);

    float sin_azimuth, cos_azimuth;
    __sincosf(atan2f(dir.y, dir.x), &sin_azimuth, &cos_azimuth);

    float sin_phi, cos_phi;
    __sincosf(2.0f*pi*curand_uniform(&rand_state), &sin_phi, &cos_phi);

    const float3 unit_v = make_float3(
        dir.z*cos_azimuth,
        dir.z*sin_azimuth,
        -sqrtf(fmaxf(0.0f, 1.0f-dir.z*dir.z))
    );
    const float3 unit_u = cross_product(unit_v, dir);

    pstruct.K_energy_dev_p[primary_idx] = K-omega;
    pstruct.dir_x_dev_p[primary_idx] = dir.x*cos_alpha+(unit_u.x*cos_phi+unit_v.x*sin_phi)*sin_alpha;
    pstruct.dir_y_dev_p[primary_idx] = dir.y*cos_alpha+(unit_u.y*cos_phi+unit_v.y*sin_phi)*sin_alpha;
    pstruct.dir_z_dev_p[primary_idx] = dir.z*cos_alpha+(unit_u.z*cos_phi+unit_v.z*sin_phi)*sin_alpha;

    pstruct.K_energy_dev_p[secondary_idx] = fermi+omega-binding;
    pstruct.dir_x_dev_p[secondary_idx] = dir.x*cos_beta-(unit_u.x*cos_phi+unit_v.x*sin_phi)*sin_beta;
    pstruct.dir_y_dev_p[secondary_idx] = dir.y*cos_beta-(unit_u.y*cos_phi+unit_v.y*sin_phi)*sin_beta;
    pstruct.dir_z_dev_p[secondary_idx] = dir.z*cos_beta-(unit_u.z*cos_phi+unit_v.z*sin_phi)*sin_beta;

    rand_state_dev_p[primary_idx] = rand_state;
}