/**
 * @file src/cdsem/octree.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED
#define eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED

#include <functional>
#include <memory>
#include <vector>
#include "point3.hh"
#include "triangle.hh"

class octree {
public:
    using triangle_p_vector = std::vector<std::shared_ptr<const triangle>>;
    octree(const point3& min, const point3& max);
    octree(const octree&);
    octree& operator=(const octree&) = delete;
    ~octree();
    std::shared_ptr<const triangle> insert(const triangle&);
    bool leaf() const;
    bool empty() const;
    int count() const;
    int depth() const;
    int occupancy() const;
    const point3& center() const;
    const point3& size() const;
    octree& root();
    const octree& root() const;
    const octree& parent() const;
    const octree& traverse(const point3&) const;
    void traverse(std::function<void(const octree&)>) const;
    triangle_p_vector::const_iterator cbegin() const;
    triangle_p_vector::const_iterator cend() const;
private:
    static const int _max_depth = 12;
    point3 _center;
    point3 _size;
    int _level;
    octree* _parent_p;
    octree* _child_p[8];
    triangle_p_vector _triangle_p_vec;
};

#endif