/**
 * @file src/cdsem/trigrid.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "trigrid.hh"
#include <algorithm>
#include "tribox.hh"

trigrid::trigrid(const point3& min, const point3& max, const index3& dim) {
    _grid_min = point3(std::min(min.x, max.x), std::min(min.y, max.y), std::min(min.z, max.z));
    _grid_max = point3(std::max(min.x, max.x), std::max(min.y, max.y), std::max(min.z, max.z));
    _grid_dim = index3(std::abs(dim.x), std::abs(dim.y), std::abs(dim.z));
}

bool trigrid::empty() const {
    return _grid_map.empty();
}

void trigrid::clear() {
    _grid_map.clear();
}

point3 trigrid::min() const {
    return _grid_min;
}

point3 trigrid::max() const {
    return _grid_max;
}

index3 trigrid::dim() const {
    return _grid_dim;
}

point3 trigrid::edge() const {
    return point3(_grid_max.x-_grid_min.x, _grid_max.y-_grid_min.y, _grid_max.z-_grid_min.z);
}

const trimesh* trigrid::at(const index3& i) const {
    if((i.x < 0) || (i.y < 0) || (i.z < 0))
        return nullptr;
    if((i.x >= _grid_dim.x) || (i.y >= _grid_dim.y) || (i.z >= _grid_dim.z))
        return nullptr;
    const int gid = i.x+_grid_dim.x*i.y+_grid_dim.x*_grid_dim.y*i.z;
    const auto cit = _grid_map.find(gid);
    if(cit != _grid_map.cend())
        return &(cit->second);
    return nullptr;
}

void trigrid::push(const trimesh& triangle_mesh, const double eps) {
    const point3 delta(edge().x/_grid_dim.x, edge().y/_grid_dim.y, edge().z/_grid_dim.z);
    for(auto cit = triangle_mesh.cbegin(); cit != triangle_mesh.cend(); cit++) {
        index3 min, max;
        std::tie(min.x, max.x) = std::minmax({
            std::floor((cit->A.x-_grid_min.x)/delta.x),
            std::floor((cit->B.x-_grid_min.x)/delta.x),
            std::floor((cit->C.x-_grid_min.x)/delta.x)
        });
        std::tie(min.y, max.y) = std::minmax({
            std::floor((cit->A.y-_grid_min.y)/delta.y),
            std::floor((cit->B.y-_grid_min.y)/delta.y),
            std::floor((cit->C.y-_grid_min.y)/delta.y)
        });
        std::tie(min.z, max.z) = std::minmax({
            std::floor((cit->A.z-_grid_min.z)/delta.z),
            std::floor((cit->B.z-_grid_min.z)/delta.z),
            std::floor((cit->C.z-_grid_min.z)/delta.z)
        });
        if((max.x < 0) || (max.y < 0) || (max.z < 0))
            continue;
        if((min.x >= _grid_dim.x) || (min.y >= _grid_dim.y) || (min.z >= _grid_dim.z))
            continue;
        for(int iz = std::max(0, min.z-1); iz <= std::min(_grid_dim.z-1, max.z+1); iz++)
        for(int iy = std::max(0, min.y-1); iy <= std::min(_grid_dim.y-1, max.y+1); iy++)
        for(int ix = std::max(0, min.x-1); ix <= std::min(_grid_dim.x-1, max.x+1); ix++) {
            double boxcenter[3];
            double boxhalfsize[3];
            double triverts[3][3];
            boxcenter[0] = _grid_min.x+(0.5+ix)*delta.x;
            boxcenter[1] = _grid_min.y+(0.5+iy)*delta.y;
            boxcenter[2] = _grid_min.z+(0.5+iz)*delta.z;
            boxhalfsize[0] = (1.0+eps)*delta.x/2;
            boxhalfsize[1] = (1.0+eps)*delta.y/2;
            boxhalfsize[2] = (1.0+eps)*delta.z/2;
            triverts[0][0] = cit->A.x; triverts[0][1] = cit->A.y; triverts[0][2] = cit->A.z;
            triverts[1][0] = cit->B.x; triverts[1][1] = cit->B.y; triverts[1][2] = cit->B.z;
            triverts[2][0] = cit->C.x; triverts[2][1] = cit->C.y; triverts[2][2] = cit->C.z;
            if(triBoxOverlap(boxcenter, boxhalfsize, triverts) > 0) {
                const int gid = ix+_grid_dim.x*iy+_grid_dim.x*_grid_dim.y*iz;
                _grid_map[gid].push(*cit);
            }
        }
    }
}

void trigrid::for_each_cell(std::function<void(const index3&, const trimesh&)> callback) const {
    for(auto cit = _grid_map.cbegin(); cit != _grid_map.cend(); cit++) {
        index3 i;
        i.x = cit->first%_grid_dim.x;
        i.y = (cit->first/_grid_dim.x)%_grid_dim.y;
        i.z = cit->first/(_grid_dim.x*_grid_dim.y);
        if(!cit->second.empty())
            callback(i, cit->second);
    }
}