/**
 * @file src/cdsem/octree.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED
#define eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED

#include <cinttypes>
#include <functional>
#include <vector>
#include "point3.hh"
#include "triangle.hh"

class octree {
public:
    using triangle_p_vector = std::vector<const::triangle*>;

    octree(const point3& min, const point3& max);
    octree(const octree&);
    octree& operator=(const octree&) = delete;
    ~octree();

    const point3& center() const;
    const point3& size() const;
    bool leaf() const;
    bool empty() const;
    int count() const;

    int level() const;
    int depth() const;
    int occupancy() const;
    uint64_t location() const;

    const triangle* insert(const triangle&);
    std::pair<double,const triangle*> intersect(const point3& A, const point3& B) const;

    octree* root();
    const octree* root() const;
    const octree* parent() const;
    const octree* traverse(int octant) const;

    triangle_p_vector::const_iterator cbegin() const;
    triangle_p_vector::const_iterator cend() const;

    int octant(const point3& pos) const;
    bool overlap(const triangle&) const;
    bool overlap(const point3& A, const point3& B) const;

public:
    static const int _max_count = 16;
    static const int _max_depth = 21;
    point3 _center;
    point3 _size;
    octree* _parent_p;
    octree* _child_p[8];
    triangle_p_vector _triangle_p_vec;
};

#endif