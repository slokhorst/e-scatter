/**
 * @file src/cdsem/octree.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#include "octree.hh"
#include <cmath>
#include <stack>
#include "point3.hh"
#include "tribox.hh"

octree::octree(const point3& min, const point3& max) {
    _center.x = (min.x+max.x)/2;
    _center.y = (min.y+max.y)/2;
    _center.z = (min.z+max.z)/2;
    _size.x = std::fabs(max.x-min.x);
    _size.y = std::fabs(max.y-min.y);
    _size.z = std::fabs(max.z-min.z);
    _level = 0;
    _parent_p = nullptr;
    for(int i = 0; i < 8; i++)
        _child_p[i] = nullptr;
}

octree::octree(const octree& node) : octree(node._center, node._size) {
    for(auto cit = node._triangle_p_vec.cbegin(); cit != node._triangle_p_vec.cend(); cit++)
        insert(*std::shared_ptr<const triangle>(*cit));
}

octree::~octree() {
    for(int i = 0; i < 8; i++)
        if(_child_p[i] != nullptr)
            delete _child_p[i];
    if(_parent_p != nullptr)
        for(int i = 0; i < 8; i++)
            if(_parent_p->_child_p[i] == this)
                _parent_p->_child_p[i] = nullptr;
}

std::shared_ptr<const triangle> octree::insert(const triangle& tri) {
    auto __overlap = [](const octree* node_p, const triangle* triangle_p)->bool {
        if((node_p == nullptr) || (triangle_p == nullptr))
            return false;
        double boxcenter[3] = {node_p->_center.x, node_p->_center.y, node_p->_center.z};
        double boxhalfsize[3] = {node_p->_size.x/2, node_p->_size.y/2, node_p->_size.z/2};
        double triverts[3][3];
        triverts[0][0] = triangle_p->A.x; triverts[0][1] = triangle_p->A.y; triverts[0][2] = triangle_p->A.z;
        triverts[1][0] = triangle_p->B.x; triverts[1][1] = triangle_p->B.y; triverts[1][2] = triangle_p->B.z;
        triverts[2][0] = triangle_p->C.x; triverts[2][1] = triangle_p->C.y; triverts[2][2] = triangle_p->C.z;
        if(triBoxOverlap(boxcenter, boxhalfsize, triverts) > 0)
            return true;
        return false;
    };
    auto __new_child = [](octree* node_p, int octant)->octree* {
        if(node_p == nullptr)
            return nullptr;
        if(node_p->_child_p[octant] != nullptr)
            return node_p->_child_p[octant];
        point3 min, max;
        min.x = node_p->_center.x-node_p->_size.x/2;
        min.y = node_p->_center.y-node_p->_size.y/2;
        min.z = node_p->_center.z-node_p->_size.z/2;
        if((octant&1) == 1)
            min.x = node_p->_center.x;
        if((octant&2) == 2)
            min.y = node_p->_center.y;
        if((octant&4) == 4)
            min.z = node_p->_center.z;
        max.x = min.x+node_p->_size.x/2;
        max.y = min.y+node_p->_size.y/2;
        max.z = min.z+node_p->_size.z/2;
        node_p->_child_p[octant] = new octree(min, max);
        node_p->_child_p[octant]->_level = node_p->_level+1;
        node_p->_child_p[octant]->_parent_p = node_p;
        return node_p->_child_p[octant];
    };
    if(!__overlap(&root(), &tri))
        return nullptr;
    std::shared_ptr<const triangle> new_triangle_p(new triangle(tri));
    std::stack<std::pair<octree*,std::shared_ptr<const triangle>>> node_p_stack;
    node_p_stack.push(std::make_pair(&root(), new_triangle_p));
    while(!node_p_stack.empty()) {
        octree* node_p = node_p_stack.top().first;
        std::shared_ptr<const triangle> triangle_p = node_p_stack.top().second;
        node_p_stack.pop();
        node_p->_triangle_p_vec.push_back(triangle_p);
        if(!node_p->leaf()) {
            /* traverse */
            for(int i = 0; i < 8; i++) {
                octree* child_p = node_p->_child_p[i];
                if(child_p == nullptr)
                    child_p = __new_child(node_p, i);
                if(__overlap(child_p, triangle_p.get()))
                    node_p_stack.push(std::make_pair(child_p, triangle_p));
            }
        } else if((node_p->count() > 8) && (node_p->_level < _max_depth)) {
            /* redistribute */
            for(int i = 0; i < 8; i++) {
                octree* child_p = node_p->_child_p[i];
                if(child_p == nullptr)
                    child_p = __new_child(node_p, i);
                for(auto cit = node_p->cbegin(); cit != node_p->cend(); cit++)
                    if(__overlap(child_p, cit->get()))
                        node_p_stack.push(std::make_pair(child_p, *cit));
            }
        }
    }
    return new_triangle_p;
}

bool octree::leaf() const {
    for(int i = 0; i < 8; i++)
        if(_child_p[i] != nullptr)
            return false;
    return true;
}

bool octree::empty() const {
    return _triangle_p_vec.empty();
}

int octree::count() const {
    return _triangle_p_vec.size();
}

int octree::depth() const {
    int _depth = 0;
    traverse([&_depth](const octree& node) {
        _depth = std::max(_depth, node._level);
    });
    return _depth;
}

int octree::occupancy() const {
    int _occupancy = 0;
    traverse([&_occupancy](const octree& node) {
        if(node.leaf())
            _occupancy = std::max(_occupancy, node.count());
    });
    return _occupancy;
}

const point3& octree::center() const {
    return _center;
}

const point3& octree::size() const {
    return _size;
}

const octree& octree::root() const {
    return root();
}

octree& octree::root() {
    octree* node_p = this;
    while(node_p->_parent_p != nullptr)
        node_p = node_p->_parent_p;
    return *node_p;
}

const octree& octree::parent() const {
    if(_parent_p == nullptr)
        return *this;
    return *_parent_p;
}

const octree& octree::traverse(const point3& pos) const {
    int i = 0;
    if(pos.x > _center.x)
        i += 1;
    if(pos.y > _center.y)
        i += 2;
    if(pos.z > _center.z)
        i += 4;
    const octree* child_p = _child_p[i];
    if(child_p == nullptr)
        return *this;
    return *child_p;
}

void octree::traverse(std::function<void(const octree&)> callback) const {
    std::stack<const octree*> node_p_stack;
    node_p_stack.push(this);
    while(!node_p_stack.empty()) {
        const octree* node_p = node_p_stack.top();
        node_p_stack.pop();
        for(int i = 0; i < 8; i++) {
            const octree* child_p = node_p->_child_p[i];
            if(child_p != nullptr)
                node_p_stack.push(child_p);
        }
        callback(*node_p);
    }
}

octree::triangle_p_vector::const_iterator octree::cbegin() const {
    return _triangle_p_vec.cbegin();
}

octree::triangle_p_vector::const_iterator octree::cend() const {
    return _triangle_p_vec.cend();
}