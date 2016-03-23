/**
 * @file src/cdsem/octree.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED
#define eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED

#include <vector>
#include "index3.hh"
#include "triangle.hh"

class octree {
public:
    struct node;
    octree(const point3& min, const point3& max);
    octree(const octree&);
    ~octree();
    octree& operator=(const octree&);
    void clear();
    int count() const;      // count the number of triangles in the tree
    int depth() const;      // determine the maximum depth of the tree
    int capacity() const;   // determine the maximum number of triangles in a leaf
    const triangle* insert(const triangle&);
    const node* traverse(const point3& pos, const node* = nullptr) const;
    std::pair<double,index3> adjacent(const node*, const point3& pos, const point3& dir) const;
    const node* adjacent(const node*, const point3& pos, const index3& dir) const;
    std::pair<double,const triangle*> intersect(const node*, const point3& pos, const point3& dir, double eps = 1e-6) const;
private:
    const int _max_depth = 10;
    std::vector<node*> _node_p_vec;
    std::vector<const triangle*> _triangle_p_vec;
};

#endif