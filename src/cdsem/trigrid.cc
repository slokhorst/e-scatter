/**
 * @file src/cdsem/trigrid.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#include "trigrid.hh"
#include <algorithm>
#include <cmath>
#include "tribox.hh"

trigrid::trigrid(const point3& min, const point3& max, const index3& dim) {
    _grid_min = point3(std::min(min.x, max.x), std::min(min.y, max.y), std::min(min.z, max.z));
    _grid_max = point3(std::max(min.x, max.x), std::max(min.y, max.y), std::max(min.z, max.z));
    _grid_dim = index3(std::abs(dim.x), std::abs(dim.y), std::abs(dim.z));
    _grid_map.assign(_grid_dim.x*_grid_dim.y*_grid_dim.z, triangle_p_vector());
}

trigrid::trigrid(const trigrid& grid) {
    *this = grid;
}

trigrid& trigrid::operator=(const trigrid& rhs) {
    if(this != &rhs) {
        clear();
        for(auto cit = rhs._grid_map.cbegin(); cit != rhs._grid_map.cend(); cit++)
            for(const triangle* triangle_p : *cit)
                insert(*triangle_p);
    }
    return *this;
}

trigrid::~trigrid() {
    clear();
}

const triangle* trigrid::insert(const triangle& tri) {
    auto __overlap = [](const point3& center, const point3& size, const triangle* triangle_p) {
        if(triangle_p == nullptr)
            return false;
        double boxcenter[3] = {center.x, center.y, center.z};
        double boxhalfsize[3] = {size.x/2, size.y/2, size.z/2};
        double triverts[3][3];
        triverts[0][0] = triangle_p->A.x; triverts[0][1] = triangle_p->A.y; triverts[0][2] = triangle_p->A.z;
        triverts[1][0] = triangle_p->B.x; triverts[1][1] = triangle_p->B.y; triverts[1][2] = triangle_p->B.z;
        triverts[2][0] = triangle_p->C.x; triverts[2][1] = triangle_p->C.y; triverts[2][2] = triangle_p->C.z;
        if(triBoxOverlap(boxcenter, boxhalfsize, triverts) > 0)
            return true;
        return false;
    };
    if(!__overlap(center(), size(), &tri))
        return nullptr;
    const triangle* triangle_p = new triangle(tri);
    _triangle_p_vec.push_back(triangle_p);
    index3 min, max;
    std::tie(min.x, max.x) = std::minmax({
        std::floor((triangle_p->A.x-_grid_min.x)/delta().x),
        std::floor((triangle_p->B.x-_grid_min.x)/delta().x),
        std::floor((triangle_p->C.x-_grid_min.x)/delta().x)
    });
    std::tie(min.y, max.y) = std::minmax({
        std::floor((triangle_p->A.y-_grid_min.y)/delta().y),
        std::floor((triangle_p->B.y-_grid_min.y)/delta().y),
        std::floor((triangle_p->C.y-_grid_min.y)/delta().y)
    });
    std::tie(min.z, max.z) = std::minmax({
        std::floor((triangle_p->A.z-_grid_min.z)/delta().z),
        std::floor((triangle_p->B.z-_grid_min.z)/delta().z),
        std::floor((triangle_p->C.z-_grid_min.z)/delta().z)
    });
    for(int iz = std::max(0, min.z-1); iz <= std::min(_grid_dim.z-1, max.z+1); iz++)
    for(int iy = std::max(0, min.y-1); iy <= std::min(_grid_dim.y-1, max.y+1); iy++)
    for(int ix = std::max(0, min.x-1); ix <= std::min(_grid_dim.x-1, max.x+1); ix++) {
        const point3 center(_grid_min.x+(0.5+ix)*delta().x, _grid_min.y+(0.5+iy)*delta().y, _grid_min.z+(0.5+iz)*delta().z);
        if(__overlap(center, delta(), triangle_p)) {
            const int gid = ix+_grid_dim.x*iy+_grid_dim.x*_grid_dim.y*iz;
            _grid_map[gid].push_back(triangle_p);
        }
    }
    return triangle_p;
}

point3 trigrid::center() const {
    return point3((_grid_min.x+_grid_max.x)/2, (_grid_min.y+_grid_max.y)/2, (_grid_min.z+_grid_max.z)/2);
}

point3 trigrid::size() const {
    return point3(_grid_max.x-_grid_min.x, _grid_max.y-_grid_min.y, _grid_max.z-_grid_min.z);
}

const index3& trigrid::dim() const {
    return _grid_dim;
}

point3 trigrid::delta() const {
    return point3(size().x/_grid_dim.x, size().y/_grid_dim.y, size().z/_grid_dim.z);
}

void trigrid::clear() {
    for(auto cit = _triangle_p_vec.cbegin(); cit != _triangle_p_vec.cend(); cit++)
        delete *cit;
    _grid_map.clear();
}
bool trigrid::empty() const {
    for(auto cit = _grid_map.cbegin(); cit != _grid_map.cend(); cit++)
        if(!cit->empty())
            return false;
    return true;
}

int trigrid::count() const {
    return _triangle_p_vec.size();
}

int trigrid::occupancy() const {
    int _occupancy = 0;
    for(auto cit = _grid_map.cbegin(); cit != _grid_map.cend(); cit++)
        _occupancy = std::max(_occupancy, static_cast<int>(cit->size()));
    return _occupancy;
}

const trigrid::triangle_p_vector& trigrid::operator[](const index3& i) const {
    const int gid = i.x+_grid_dim.x*i.y+_grid_dim.x*_grid_dim.y*i.z;
    return _grid_map.at(gid);
}

void trigrid::for_each_cell(std::function<void(const index3&, const triangle_p_vector&)> callback) const {
    for(size_t i = 0; i < _grid_map.size(); i++)
        if(!_grid_map[i].empty()) {
            index3 j;
            j.x = i%_grid_dim.x;
            j.y = (i/_grid_dim.x)%_grid_dim.y;
            j.z = (i/(_grid_dim.x*_grid_dim.y))%_grid_dim.z;
            callback(j, _grid_map[i]);
        }
}