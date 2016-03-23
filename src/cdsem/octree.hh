/**
 * @file src/cdsem/octree.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED
#define eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED

#include <vector>
#include "triangle.hh"

class octree {
public:
    octree(const point3& min_box, const point3& max_box);
    octree(const octree&);
    ~octree();
    octree& operator=(const octree&);
    void clear();
    int count() const;      // count the number of triangles in the tree
    int depth() const;      // determine the maximum depth of the tree
    int capacity() const;   // determine the maximum number of triangles in a leaf
    const triangle* insert(const triangle&);
private:
    const int _max_depth = 10;
    struct node;
    std::vector<node*> _node_p_vec;
    std::vector<const triangle*> _triangle_p_vec;
};

#endif