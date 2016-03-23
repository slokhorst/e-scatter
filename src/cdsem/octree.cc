/**
 * @file src/cdsem/octree.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#include "octree.hh"
#include <algorithm>
#include <cmath>
#include <set>
#include <stack>
#include "point3.hh"
#include "tribox.hh"

struct octree::node {
    bool overlap(const triangle& _triangle) const {
        double boxcenter[3] = {center.x, center.y, center.z};
        double boxhalfsize[3] = {size.x/2, size.y/2, size.z/2};
        double triverts[3][3];
        triverts[0][0] = _triangle.A.x; triverts[0][1] = _triangle.A.y; triverts[0][2] = _triangle.A.z;
        triverts[1][0] = _triangle.B.x; triverts[1][1] = _triangle.B.y; triverts[1][2] = _triangle.B.z;
        triverts[2][0] = _triangle.C.x; triverts[2][1] = _triangle.C.y; triverts[2][2] = _triangle.C.z;
        if(triBoxOverlap(boxcenter, boxhalfsize, triverts) > 0)
            return true;
        return false;
    }
    bool leaf() const {
        for(int i = 0; i < 8; i++)
            if(child_p[i] != nullptr)
                return false;
        return true;
    }
    point3 center;
    point3 size;
    int depth;
    node* parent_p = nullptr;
    node* child_p[8] = {nullptr};
    std::vector<const triangle*> triangle_p_vec;
};

octree::octree(const point3& min, const point3& max) {
    node* root_p = new node();
    _node_p_vec.push_back(root_p);
    root_p->center.x = (min.x+max.x)/2;
    root_p->center.y = (min.y+max.y)/2;
    root_p->center.z = (min.z+max.z)/2;
    root_p->size.x = std::fabs(max.x-min.x);
    root_p->size.y = std::fabs(max.y-min.y);
    root_p->size.z = std::fabs(max.z-min.z);
    root_p->depth = 0;
}

octree::octree(const octree& tree) {
    *this = tree;
}

octree& octree::operator=(const octree& rhs) {
    if(this != &rhs) {
        clear();
        node* root_p = new node();
        _node_p_vec.push_back(root_p);
        root_p->center = rhs._node_p_vec.front()->center;
        root_p->size = rhs._node_p_vec.front()->size;
        root_p->depth = 0;
        for(auto cit = rhs._triangle_p_vec.cbegin(); cit != rhs._triangle_p_vec.cend(); cit++) {
            const triangle* triangle_p = *cit;
            insert(*triangle_p);
        }
    }
    return *this;
}

octree::~octree() {
    clear();
}

void octree::clear() {
    for(auto it = _node_p_vec.begin(); it != _node_p_vec.end(); it++)
        delete *it;
    _node_p_vec.clear();
    for(auto it = _triangle_p_vec.begin(); it != _triangle_p_vec.end(); it++)
        delete *it;
    _triangle_p_vec.clear();
}

int octree::count() const {
    const node* root_p = _node_p_vec.front();
    std::set<const triangle*> triangle_p_set;
    std::stack<const node*> node_p_stack;
    node_p_stack.push(root_p);
    while(!node_p_stack.empty()) {
        const node* node_p = node_p_stack.top();
        node_p_stack.pop();
        for(auto cit = node_p->triangle_p_vec.cbegin(); cit != node_p->triangle_p_vec.cend(); cit++)
            triangle_p_set.insert(*cit);
        for(int i = 0; i < 8; i++) {
            const node* child_p = node_p->child_p[i];
            if(child_p != nullptr)
                node_p_stack.push(child_p);
        }
    }
    return triangle_p_set.size();
}

int octree::depth() const {
    int max_depth = 0;
    const node* root_p = _node_p_vec.front();
    std::stack<const node*> node_p_stack;
    node_p_stack.push(root_p);
    while(!node_p_stack.empty()) {
        const node* node_p = node_p_stack.top();
        node_p_stack.pop();
        max_depth = std::max(max_depth, node_p->depth);
        for(int i = 0; i < 8; i++) {
            const node* child_p = node_p->child_p[i];
            if(child_p != nullptr)
                node_p_stack.push(child_p);
        }
    }
    return max_depth;
}

int octree::capacity() const {
    int max_capacity = 0;
    const node* root_p = _node_p_vec.front();
    std::stack<const node*> node_p_stack;
    node_p_stack.push(root_p);
    while(!node_p_stack.empty()) {
        const node* node_p = node_p_stack.top();
        node_p_stack.pop();
        max_capacity = std::max(max_capacity, static_cast<int>(node_p->triangle_p_vec.size()));
        for(int i = 0; i < 8; i++) {
            const node* child_p = node_p->child_p[i];
            if(child_p != nullptr)
                node_p_stack.push(child_p);
        }
    }
    return max_capacity;
}

const triangle* octree::insert(const triangle& tri) {
    node* root_p = _node_p_vec.front();
    if(!root_p->overlap(tri))
        return nullptr;
    const triangle* triangle_p = new triangle(tri);
    _triangle_p_vec.push_back(triangle_p);
    std::stack<std::pair<node*,const triangle*>> node_p_stack;
    node_p_stack.push(std::make_pair(root_p, triangle_p));
    while(!node_p_stack.empty()) {
        node* node_p = node_p_stack.top().first;
        triangle_p = node_p_stack.top().second;
        node_p_stack.pop();
        if(!node_p->leaf()) {
            /* traverse */
            for(int i = 0; i < 8; i++) {
                node* child_p = node_p->child_p[i];
                if(child_p->overlap(*triangle_p))
                    node_p_stack.push(std::make_pair(child_p, triangle_p));
            }
        } else {
            /* push triangle */
            node_p->triangle_p_vec.push_back(triangle_p);
            if((node_p->triangle_p_vec.size() > 8) && (node_p->depth < _max_depth)) {
                /* redistribute triangles */
                for(int i = 0; i < 8; i++) {
                    node* child_p = new node();
                    _node_p_vec.push_back(child_p);
                    node_p->child_p[i] = child_p;
                    child_p->center.x = node_p->center.x-node_p->size.x/4;
                    child_p->center.y = node_p->center.y-node_p->size.y/4;
                    child_p->center.z = node_p->center.z-node_p->size.z/4;
                    if((i&1) == 1)
                        child_p->center.x = node_p->center.x+node_p->size.x/4;
                    if((i&2) == 2)
                        child_p->center.y = node_p->center.y+node_p->size.y/4;
                    if((i&4) == 4)
                        child_p->center.z = node_p->center.z+node_p->size.z/4;
                    child_p->size.x = node_p->size.x/2;
                    child_p->size.y = node_p->size.y/2;
                    child_p->size.z = node_p->size.z/2;
                    child_p->depth = node_p->depth+1;
                    child_p->parent_p = node_p;
                    for(auto cit = node_p->triangle_p_vec.cbegin(); cit != node_p->triangle_p_vec.cend(); cit++)
                        if(child_p->overlap(**cit))
                            node_p_stack.push(std::make_pair(child_p, *cit));
                }
                node_p->triangle_p_vec.clear();
            }
        }
    }
    return _triangle_p_vec.back();
}