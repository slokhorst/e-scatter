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

/*! \brief Compute unit vector for Euler angles
 *  @return
 *    Float3 of length 1
 */
__device__ float3 make_unit_vec(float cos_theta, float phi) {
    float sin_phi, cos_phi;
    sincosf(phi, &sin_phi, &cos_phi);
    cos_theta = clamp(cos_theta, -1.0f, 1.0f);
    const float sin_theta = sqrtf(1.0f-cos_theta*cos_theta);
    return make_float3(sin_theta*cos_phi, sin_theta*sin_phi, cos_theta);
}

/*! \brief Determination of a normal vector
 *   This function is only used to determine a random scatter angle for events.
 *   The angle phi is with reference to a fixed non-determined vector.
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

/*! \brief Point inside AABB (Axis Aligned Bounding Box) test.
 *  @param[in] pos
 *  @param[in] center
 *  @param[in] halfsize
 *  @return
 *    True if point is inside, but not on the boundary.
 */
__device__ bool inside_AABB(float3 pos, float3 center, float3 halfsize) {
    if((pos.x > center.x-halfsize.x) && (pos.x < center.x+halfsize.x))
    if((pos.y > center.y-halfsize.y) && (pos.y < center.y+halfsize.y))
    if((pos.z > center.z-halfsize.z) && (pos.z < center.z+halfsize.z))
        return true;
    return false;
}

/*! \brief Determine time of interesection with each orthogonal plane of the AABB.
 *   intersect.{x,y,z} = pos.{x,y,z} + time.{x,y,z}*dir.{x,y,z}
 *   Where time.{x,y,z} is the corresponding component of the returned Float3.
 *  @param[in] pos
 *  @param[in] dir
 *  @param[in] center
 *  @param[in] halfsize
 *  @return
 *    Float3 of relative orthogonal plane intersections.
 */
__device__ float3 AABB_intersect(float3 pos, float3 dir, float3 center, float3 halfsize) {
    return make_float3(
        (center.x+copysignf(halfsize.x+_eps, dir.x)-pos.x)/dir.x,
        (center.y+copysignf(halfsize.y+_eps, dir.y)-pos.y)/dir.y,
        (center.z+copysignf(halfsize.z+_eps, dir.z)-pos.z)/dir.z
    );
}

/*! \brief Determine ray-triangle interesection.
 *   See: T. Möller and B. Trumbore, Journal of Graphics Tools, 2(1):21--28, 1997.
 *  @param[in] pos
 *  @param[in] dir
 *  @param[in] O point of origin for the triangle.
 *  @param[in] e1 vector of first edge
 *  @param[in] e2 vector of second edge
 *  @note
      Given a triangle with vertices A, B and C. The origin is A,
      and the first edge is B-A and the second edge is C-A. The algorithm
      is insensitive to orientation of the triangle.
 *  @return
 *    Float respresenting the time of intersection OR -1.0f if no intersection.
 */
__device__ float triangle_intersect(float3 pos, float3 dir, float3 O, float3 e1, float3 e2) {
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

/*! \brief Wrapper for curand uniform random number generator.
 *    Reads a random state from global memory, allows to sample random values and
 *    upon destruction updates the random state in global memory.
 */
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

/*! \brief CUDA kernel to initialize random states in global memory.
 *   For more detailed information, see http://docs.nvidia.com/cuda/curand/
 *  @note The seed for each individial thread is determined by the global seed
 *        and the thread index. The default random number generator used in curand
 *        is the XORWOW generator.
 */
__global__ void cuda_init_rand_state(curandState* rand_state_p, unsigned long long seed, int n, scatter_options opt) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= n)
        return;
    curand_init(seed, i, 0, &(rand_state_p[i]));
}

/*! \brief Determine the distance (in nm) to travel for a particular event.
 *  writes to members of pstruct: status_dev_p, distance_dev_p and rand_state_dev_p
 */
__global__ void cuda_init_trajectory(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p, scatter_options opt) {
    const int particle_idx = threadIdx.x+blockIdx.x*blockDim.x;
    if(particle_idx >= pstruct.capacity)
        return;

    // ignore all particles that are TERMINATED, DETECTED or PENDING.
    switch(pstruct.status_dev_p[particle_idx]) {
        case cuda_particle_struct::TERMINATED:
        case cuda_particle_struct::DETECTED:
        case cuda_particle_struct::PENDING:
            return;
        default:
            break;
    }

    // retrieve particle position from global memory
    const float3 pos = make_float3(
        pstruct.pos_x_dev_p[particle_idx],
        pstruct.pos_y_dev_p[particle_idx],
        pstruct.pos_z_dev_p[particle_idx]
    );

    // KILL and EXIT if particle is not inside the root cell of the octree
    if(!inside_AABB(pos, gstruct.AABB_center, gstruct.AABB_halfsize)) {
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
        return;
    }

    // retrieve the current material index of the particle
    // note that a negative material index implies vacuum
    const int material_idx = pstruct.material_idx_dev_p[particle_idx];
    if(material_idx < 0) {
        // vacuum propgation
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::NO_EVENT;
        pstruct.distance_dev_p[particle_idx] =
            2.0f*norm3df(gstruct.AABB_halfsize.x, gstruct.AABB_halfsize.y, gstruct.AABB_halfsize.z);
        return;
    }
    // ELSE: material propagation...

    // retrieve the kinetic energy (in eV) of the particle from global memory
    const float K = pstruct.K_energy_dev_p[particle_idx];

    // KILL and EXIT: if energy is below potential barrier
    if(K < mstruct.barrier_dev_p[material_idx]) {
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
        return;
    }
    // ELSE: particle may escape the material...

    // determine the total inverse mean free path by log-log interpolation of
    // the elastic and inelastic inverse mean free path
    // TODO: put in inline function
    float elastic_imfp, inelastic_imfp, total_imfp;
    {
        const float x = logf(K/mstruct.K_min)/logf(mstruct.K_max/mstruct.K_min);
        const int3 offset = make_int3(0, 0, material_idx);
        const int dim = mstruct.K_cnt;
        elastic_imfp = expf(interp1(mstruct.elastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, offset, dim, x));
        inelastic_imfp = expf(interp1(mstruct.inelastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, offset, dim, x));
        total_imfp = elastic_imfp+inelastic_imfp;
    }

    curand_wrapper curand(rand_state_dev_p[particle_idx]);

    // draw a random distance and fix the event (see thesis T.V. Eq. 3.8)
    pstruct.distance_dev_p[particle_idx] = -logf(curand.uniform())/total_imfp;
    if(curand.uniform() < elastic_imfp/total_imfp)
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::ELASTIC_EVENT;
    else
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::INELASTIC_EVENT;
}

/*! \brief Attempts to displace the particle a distance `distance_dev_p`.
 *         If a triangle intersection occurs, the `status_dev_p` is updated accordingly.
 *         In this case the triangle index `triangle_idx_dev_p` is set to blocking triangle.
 *         The particle position is updated to the new position, which is one of the following three:
 *           (1) pos+dir*distance in the case of no intersection.
 *           (2) point of intersection with the blocking triangle.
 *           (3) orthogonal plane of the root cell -> particle is terminated.
 */
__global__ void cuda_update_trajectory(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, scatter_options opt) {
    const int particle_idx = threadIdx.x+blockIdx.x*blockDim.x;
    if(particle_idx >= pstruct.capacity)
        return;

    // ignore all particles that are TERMINATED, DETECTED or PENDING.
    switch(pstruct.status_dev_p[particle_idx]) {
        case cuda_particle_struct::TERMINATED:
        case cuda_particle_struct::DETECTED:
        case cuda_particle_struct::PENDING:
            return;
        default:
            break;
    }

    // retrieve particle position (in nm) from global memory
    float3 pos = make_float3(
        pstruct.pos_x_dev_p[particle_idx],
        pstruct.pos_y_dev_p[particle_idx],
        pstruct.pos_z_dev_p[particle_idx]
    );

    // retrieve particle direction from global memory and normalize
    float3 dir = make_float3(
        pstruct.dir_x_dev_p[particle_idx],
        pstruct.dir_y_dev_p[particle_idx],
        pstruct.dir_z_dev_p[particle_idx]
    );
    dir *= rnorm3df(dir.x, dir.y, dir.z);

    // retrieve current material index (guaranteed to be larger than 0)
    const int material_idx = pstruct.material_idx_dev_p[particle_idx];

    // retrieve the last triangle of interface scattering from global memory
    const int last_triangle_idx = pstruct.triangle_idx_dev_p[particle_idx];

    // retrieve the kinetic energy (in eV) from global memory
    const float K = pstruct.K_energy_dev_p[particle_idx];

    // retrieve the intended distance (in nm) to travel (as determined by `init_trajectory`)
    float distance = pstruct.distance_dev_p[particle_idx];

    // Search in octree
    // A cell in the octree is uniquely identified by a location code.
    // The root cell has location code 1.
    // The location code of a decedant is obtained by shifting the location code
    // to the left by three bits and adding the octant of the decadant to the location code.
    uint64_t location = 1;
    do {
        // define the root axis aligned bounding box
        float4 AABB = make_float4(
            gstruct.AABB_center.x,
            gstruct.AABB_center.y,
            gstruct.AABB_center.z,
            1.0f // scale factor
        );

        // initial index to linearized octree
        int index = 0;

        // traverse to location: we compute the location of the current node and
        //   set the index to the node in the linearized octree array
        // `__clzll` computes number of leading zeros, effectively this is 64-floor(log2(n))
        for(int i = 60-__clzll(location); i >= 0; i -= 3) {
            const unsigned int octant = (location>>i)&7;
            index = gstruct.octree_dev_p[index+octant];
            AABB.x += AABB.w*gstruct.AABB_halfsize.x*(1.0f*(octant&1)-0.5f);
            AABB.y += AABB.w*gstruct.AABB_halfsize.y*(1.0f*((octant&2)>>1)-0.5f);
            AABB.z += AABB.w*gstruct.AABB_halfsize.z*(1.0f*((octant&4)>>2)-0.5f);
            AABB.w *= 0.5f;
        }

        // traverse to leaf node (which has a linearized octree index smaller than zero)
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
        // note that a leaf has a guaranteed negative index, the actual index to the leaf
        // node is then obtained by negating.
        index = -index;

        // determine cell/triangle intersections for the leaf
        float intersect;
        int target; // target >= 0 : triangle index
                    // target = -1 : cell x intersection
                    // target = -2 : cell y intersection
                    // target = -4 : cell z intersection

        // find intersection times with each orthogonal plane of the leaf's bounding box,
        // then see which plane is reached first.
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

        // iterate triangles in leaf
        while(true) {
            const int triangle_idx = gstruct.octree_dev_p[index++];
            if(triangle_idx < 0)
                break;
            // don't intersect with last triangle
            if(triangle_idx == last_triangle_idx)
                continue;

            // retrieve triangle geometry from global memory
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

            // see on which side of triangle we are, and find the correct
            // material on its other side
            int mat_idx_out;
            if(dot_product(cross_product(e1, e2), dir) < 0)
                mat_idx_out = gstruct.material_idx_in_dev_p[triangle_idx];
            else
                mat_idx_out = gstruct.material_idx_out_dev_p[triangle_idx];

            // if the outgoing material is the same as current, nothing happens
            // if the triangle represents a detector which can't see the current
            // particle, nothing happens
            if((mat_idx_out == material_idx) || (mat_idx_out == triangle::NOP))
                continue;
            else if((mat_idx_out == triangle::DETECTOR_LT50) && (K >= 50))
                continue;
            else if((mat_idx_out == triangle::DETECTOR_GE50) && (K < 50))
                continue;

            // else, compute the intersection with the triangle; keep it if it
            // is closer than the current distance of scattering.
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
            // EXIT: triangle intersection updates status and position of the particle
            pstruct.status_dev_p[particle_idx] = cuda_particle_struct::INTERSECT_EVENT;
            pstruct.triangle_idx_dev_p[particle_idx] = target;
            pstruct.pos_x_dev_p[particle_idx] = pos.x+dir.x*intersect;
            pstruct.pos_y_dev_p[particle_idx] = pos.y+dir.y*intersect;
            pstruct.pos_z_dev_p[particle_idx] = pos.z+dir.z*intersect;
            return;
        }

        // go to edge of the leaf
        distance -= intersect;
        pos.x += dir.x*intersect;
        pos.y += dir.y*intersect;
        pos.z += dir.z*intersect;

        // find adjacent node by bit-magic, and loop
        // the adjacent node is determined solely by location code!
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

    // EXIT: particle is out of grid (left root cell)
    pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
    return;
}

// ************************************************************************** //
// remaining kernels manage the individual event types:
//   INTERSECT_EVENT -> quantum mechanical transmission/reflection/refraction (computed by analytical formulas).
//   ELASTIC_EVENT   -> ELSEPA tables and phonon scattering (from single table).
//   INELASTIC_EVENT -> Dielectric function model (Ashley with Kieft and Bosch refinements).

__global__ void cuda_intersection_event(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p, scatter_options opt) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= pstruct.capacity)
        return;
    const int particle_idx = pstruct.particle_idx_dev_p[i];

    // ignore all particles except those with an intersection event.
    switch(pstruct.status_dev_p[particle_idx]) {
        case cuda_particle_struct::INTERSECT_EVENT:
            break;
        default:
            return;
    }

    // retrieve direction from global memory and normalize
    // TODO: make inline function?
    float3 dir = make_float3(
        pstruct.dir_x_dev_p[particle_idx],
        pstruct.dir_y_dev_p[particle_idx],
        pstruct.dir_z_dev_p[particle_idx]
    );
    dir *= rnorm3df(dir.x, dir.y, dir.z);

    // retrieve triangle index and the edges of the triangle.
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

    // determine normal vector of the triangle and normalize to unity
    float3 normal = cross_product(e1, e2);
    normal *= rnorm3df(normal.x, normal.y, normal.z);

    // note __saturatef(x) equals clampf(x, 0, 1)
    const float cos_theta = dot_product(normal, dir);
    const float sin_theta = sqrtf(1.0f-__saturatef(cos_theta*cos_theta));

    // determine the material index in the following way:
    //   material_idx_in represents the current material
    //   material_idx_out represents the material when passing through the interface
    int material_idx_in, material_idx_out;
    if(cos_theta > 0) {
        material_idx_in = gstruct.material_idx_in_dev_p[triangle_idx];
        material_idx_out = gstruct.material_idx_out_dev_p[triangle_idx];
    } else {
        material_idx_in = gstruct.material_idx_out_dev_p[triangle_idx];
        material_idx_out = gstruct.material_idx_in_dev_p[triangle_idx];
    }

    // manage special cases for electron detection, electron mirrors and terminators.
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

    // determine the change in energy `dU` (in eV) when passing through the interface
    // see thesis T.V. Eq. 3.136
    float dU = 0;
    if(material_idx_in >= 0) {
        dU -= mstruct.barrier_dev_p[material_idx_in];
    }
    if(material_idx_out >= 0) {
        dU += mstruct.barrier_dev_p[material_idx_out];
    }

    curand_wrapper curand(rand_state_dev_p[particle_idx]);

    const float K = pstruct.K_energy_dev_p[particle_idx];

    // determine transmission probability (only if energy suffices)
    // see thesis T.V. Eq. 3.145
    if(K*cos_theta*cos_theta+dU > 0.0f) {
        const float s = sqrtf(1.0f+dU/(K*cos_theta*cos_theta));
        const float T = (opt.quantum_transmission_flag) ? 4.0f*s/((1.0f+s)*(1.0f+s)) : 1.0f;
        if(curand.uniform() < T) {
            if(opt.interface_refraction_flag) {
                // determine the angle of refraction
                // see thesis T.V. Eq. 3.139
                pstruct.dir_x_dev_p[particle_idx] = (dir.x-normal.x*cos_theta)+s*normal.x*cos_theta;
                pstruct.dir_y_dev_p[particle_idx] = (dir.y-normal.y*cos_theta)+s*normal.y*cos_theta;
                pstruct.dir_z_dev_p[particle_idx] = (dir.z-normal.z*cos_theta)+s*normal.z*cos_theta;
            }
            // if there is transmission, then adjust the kinetic energy,
            // update the current material index and EXIT.
            pstruct.K_energy_dev_p[particle_idx] = K+dU;
            pstruct.material_idx_dev_p[particle_idx] = material_idx_out;
            return;
        }
    }

    // surface absorption? (this is in accordance with Kieft & Bosch code)
    // note that the default behaviour has this feature disabled
    if((opt.interface_absorption_flag) && (dU < 0) && (curand.uniform() < __expf(1.0f+0.5f*K/dU))) {
        pstruct.status_dev_p[particle_idx] = cuda_particle_struct::TERMINATED;
        return;
    }

    // the only remaining case is total internal reflection
    pstruct.dir_x_dev_p[particle_idx] = dir.x-2.0f*normal.x*cos_theta;
    pstruct.dir_y_dev_p[particle_idx] = dir.y-2.0f*normal.y*cos_theta;
    pstruct.dir_z_dev_p[particle_idx] = dir.z-2.0f*normal.z*cos_theta;
}

__global__ void cuda_elastic_event(cuda_particle_struct pstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p, scatter_options opt) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= pstruct.capacity)
        return;
    const int particle_idx = pstruct.particle_idx_dev_p[i];

    // ignore all particles except those with an elastic event.
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
    {// draw a random elastic scatter angle by interpolating tables
        const float x = logf(K/mstruct.K_min)/logf(mstruct.K_max/mstruct.K_min);
        const float y = curand.uniform();
        const int3 offset = make_int3(0, 1, material_idx);
        const int2 dim = make_int2(mstruct.K_cnt, mstruct.P_cnt);
        cos_theta = interp2(mstruct.elastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, offset, dim, x, y);
        cos_theta = clamp(cos_theta, -1.0f, 1.0f);
        sin_theta = sqrtf(1.0f-cos_theta*cos_theta);
    }

    // retrieve the current direction from global memory and normalize
    float3 dir = make_float3(
        pstruct.dir_x_dev_p[particle_idx],
        pstruct.dir_y_dev_p[particle_idx],
        pstruct.dir_z_dev_p[particle_idx]
    );
    dir *= rnorm3df(dir.x, dir.y, dir.z);

    // find a random normal vector to the current direction of flight and normalize
    float3 normal_dir = make_normal_vec(dir, 2.0f*_pi*curand.uniform());
    normal_dir *= rnorm3df(normal_dir.x, normal_dir.y, normal_dir.z);

    // determine the scattered direction
    pstruct.dir_x_dev_p[particle_idx] = dir.x*cos_theta+normal_dir.x*sin_theta;
    pstruct.dir_y_dev_p[particle_idx] = dir.y*cos_theta+normal_dir.y*sin_theta;
    pstruct.dir_z_dev_p[particle_idx] = dir.z*cos_theta+normal_dir.z*sin_theta;

    // special cases for phonon scattering and atom recoil energy loss.
    // the energy domain for phonon scattering is exactly the same as in the
    //   original Kieft & Bosch code.
    // the amount of energy loss for phonons can be found in the thesis of T.V. Eq. 3.116.
    // TODO: set variable domain for phonon scattering
    if(K < 200) {
        if(opt.acoustic_phonon_loss_flag)
            pstruct.K_energy_dev_p[particle_idx] = K-fminf(50e-3f, mstruct.phonon_loss_dev_p[material_idx]);
    } else {
        // account for atomic recoil (only added for compliance with Kieft & Bosch code)
        // There is no reference for this formula, can only be found in the Kieft & Bosch code.
        // Default behaviour does not include this effect.
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

    // Ignore all particles except those with an inelastic event and those that are pending.
    // Pending particles are actually particles with an inelastic event but with
    //   no room for creating a secondary in the array.
    switch(pstruct.status_dev_p[primary_idx]) {
        case cuda_particle_struct::INELASTIC_EVENT:
        case cuda_particle_struct::PENDING:
            break;
        default:
            return;
    }

    // determine a free slot for a new particle.
    // note that this index is obtained by sorting the particle array in such a way
    //   that terminated particles reside at the end of the array.
    const int secondary_idx = pstruct.particle_idx_dev_p[pstruct.capacity-1-i];

    // change the status to pending if the particle at the back is not terminated and EXIT.
    if(pstruct.status_dev_p[secondary_idx] != cuda_particle_struct::TERMINATED) {
        pstruct.status_dev_p[primary_idx] = cuda_particle_struct::PENDING;
        return;
    }
    // otherwise: this is an inelastic event where a secondary can be created.
    pstruct.status_dev_p[primary_idx] = cuda_particle_struct::INELASTIC_EVENT;

    // forget last intersected triangle
    pstruct.triangle_idx_dev_p[primary_idx] = -1;

    // retrieve material index from global memory
    const int material_idx = pstruct.material_idx_dev_p[primary_idx];

    // retrieve the fermi energy (in eV) of the current material
    const float fermi = mstruct.fermi_dev_p[material_idx];

    // retrieve the kinetic energy (in eV) of the particle
    const float K = pstruct.K_energy_dev_p[primary_idx];

    curand_wrapper curand(rand_state_dev_p[primary_idx]);

    // draw a random zero-momentum energy loss of the primary electron
    float omega0;
    {// see thesis T.V. Eq. 3.82.
        const float x = logf(K/mstruct.K_min)/logf(mstruct.K_max/mstruct.K_min);
        const float y = curand.uniform();
        const int3 offset = make_int3(0, 1, material_idx);
        const int2 dim = make_int2(mstruct.K_cnt, mstruct.P_cnt);
        omega0 = expf(interp2(mstruct.inelastic_dev_p, mstruct.pitch, mstruct.P_cnt+1, offset, dim, x, y));
        omega0 = posf(omega0);
    }

    // draw a random binding energy of the secondary electron.
    float binding;
    {// see thesis T.V. Fig. 3.14 for an example of the possible binding energies.
        const float x = logf(omega0/mstruct.K_min)/logf(mstruct.K_max/mstruct.K_min);
        const float y = curand.uniform();
        const int ix = clamp(__float2int_rd(x*(mstruct.K_cnt-1)), 0, mstruct.K_cnt-1);
        const int iy = clamp(__float2int_rd(y*(mstruct.P_cnt-1)), 0, mstruct.P_cnt-1);
        binding = cuda_make_ptr(mstruct.ionization_dev_p, mstruct.pitch, mstruct.P_cnt+1, iy, material_idx)[ix];
    }

    // draw a random total energy loss for the primary electron
    float omega;
    {// see thesis T.V. Eq. 3.85.
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

    // special cases if there is no binding energy:
    if(binding < 0) {
        const float band_gap = mstruct.band_gap_dev_p[material_idx];
        if(band_gap < 0) {
            // TODO for metals: excitation of a fermi sea electron
        } else if(omega0 > band_gap) {
            // electron excitation across the band gap (see page 78 thesis T.V.)
            binding = band_gap;
        } else {
            // sub-band gap energy loss in semiconductors and insulators (see page 78 thesis T.V.)
            // energy loss due to longitudinal optical phonon excitation is assumed
            // update energy and EXIT
            if(opt.optical_phonon_loss_flag)
                pstruct.K_energy_dev_p[primary_idx] = K-omega0;
            return;
        }
    }
    binding = posf(binding);

    // retrieve direction from global memory and normalize
    float3 primary_dir = make_float3(
        pstruct.dir_x_dev_p[primary_idx],
        pstruct.dir_y_dev_p[primary_idx],
        pstruct.dir_z_dev_p[primary_idx]
    );
    primary_dir *= rnorm3df(primary_dir.x, primary_dir.y, primary_dir.z);

    // determine random normal vector to determine the scattering direction
    float3 normal_dir = make_normal_vec(primary_dir, 2.0f*_pi*curand.uniform());
    normal_dir *= rnorm3df(normal_dir.x, normal_dir.y, normal_dir.z);

    // Determine the inelastic scattering angles.
    // We have strictly followed the method of Ivanchenko (and thus Kieft and Bosch)
    // See for more details thesis T.V. page 80-85.
    const float _K = K-fermi+2.0f*binding; // see thesis T.V. Eq. 3.105
    const float dK = binding+omega;        // see thesis T.V. Eq. 3.106
    const float cos_theta = sqrtf(dK/_K);  // see thesis T.V. Eq. 3.100
    const float sin_theta = sqrtf(1.0f-__saturatef(cos_theta*cos_theta));

    // Determine initial secondary direction (see thesis T.V. Eq. 3.107)
    // The initial direction is determined by assuming that the secondary electron
    //   is at rest
    float3 secondary_dir = make_float3(
        primary_dir.x*cos_theta+normal_dir.x*sin_theta,
        primary_dir.y*cos_theta+normal_dir.y*sin_theta,
        primary_dir.z*cos_theta+normal_dir.z*sin_theta
    );

    // Add (optional) direction to account for the (intrinsic) instantaneous momentum
    //  of the secondary electron.
    // See thesis T.V. Eq. 3.108
    if(opt.instantaneous_momentum_flag) {
        const float random_cos = 2.0f*curand.uniform()-1.0f;
        const float random_phi = 2.0f*_pi*curand.uniform();
        secondary_dir *= rnorm3df(secondary_dir.x, secondary_dir.y, secondary_dir.z);
        secondary_dir += sqrtf(binding/dK)*make_unit_vec(random_cos, random_phi);
    }

    // ensure proper normalization of the secondary directional vector.
    secondary_dir *= rnorm3df(secondary_dir.x, secondary_dir.y, secondary_dir.z);

    if(opt.generate_secondary_flag) {
        // new secondaries receive the status of `NEW_SECONDARY`.
        pstruct.status_dev_p[secondary_idx] = cuda_particle_struct::NEW_SECONDARY;

        // the particle tag and material index is inherited from the primary electron.
        pstruct.material_idx_dev_p[secondary_idx] = material_idx;
        pstruct.particle_tag_dev_p[secondary_idx] = pstruct.particle_tag_dev_p[primary_idx];

        // the position of the secondary is also inherited from the primary electron.
        pstruct.pos_x_dev_p[secondary_idx] = pstruct.pos_x_dev_p[primary_idx];
        pstruct.pos_y_dev_p[secondary_idx] = pstruct.pos_y_dev_p[primary_idx];
        pstruct.pos_z_dev_p[secondary_idx] = pstruct.pos_z_dev_p[primary_idx];

        // determine the kinetic energy of the secondary electron
        pstruct.K_energy_dev_p[secondary_idx] = fermi+omega-binding; // See thesis T.V. Eq. 3.86
        pstruct.dir_x_dev_p[secondary_idx] = secondary_dir.x;
        pstruct.dir_y_dev_p[secondary_idx] = secondary_dir.y;
        pstruct.dir_z_dev_p[secondary_idx] = secondary_dir.z;
    }

    // Update the kinetic energy of the primary electron.
    // Note that energy is conserved.
    pstruct.K_energy_dev_p[primary_idx] = K-omega;

    // primary direction determined by non-relativistic momentum-conservation, i.e.:
    //   sin(theta)*primary_dir_2 = primary_dir - cos(theta)*secondary_dir;
    // See thesis T.V. Eq. 3.111
    if(opt.momentum_conservation_flag) {
        pstruct.dir_x_dev_p[primary_idx] = primary_dir.x-cos_theta*secondary_dir.x;
        pstruct.dir_y_dev_p[primary_idx] = primary_dir.y-cos_theta*secondary_dir.y;
        pstruct.dir_z_dev_p[primary_idx] = primary_dir.z-cos_theta*secondary_dir.z;
    }
}