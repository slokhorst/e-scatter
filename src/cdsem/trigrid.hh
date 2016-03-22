/**
 * @file src/cdsem/trigrid.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__TRIGRID__HEADER_INCLUDED
#define eSCATTER__CDSEM__TRIGRID__HEADER_INCLUDED

#include <cfloat>
#include <functional>
#include <map>
#include "index3.hh"
#include "point3.hh"
#include "trimesh.hh"

class trigrid {
public:
    trigrid(const point3& min, const point3& max, const index3& dim);
    bool empty() const;
    void clear();
    point3 min() const;
    point3 max() const;
    index3 dim() const;
    point3 edge() const;
    const trimesh* at(const index3&) const;
    void push(const trimesh&, const double eps = FLT_EPSILON);
    void for_each_cell(std::function<void(const index3&, const trimesh&)>) const;
private:
    point3 _grid_min;
    point3 _grid_max;
    index3 _grid_dim;
    std::map<int,trimesh> _grid_map;
};

#endif