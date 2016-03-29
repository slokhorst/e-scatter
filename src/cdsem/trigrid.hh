/**
 * @file src/cdsem/trigrid.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__TRIGRID__HEADER_INCLUDED
#define eSCATTER__CDSEM__TRIGRID__HEADER_INCLUDED

#include <functional>
#include <vector>
#include "index3.hh"
#include "point3.hh"
#include "triangle.hh"

class trigrid {
public:
    using triangle_p_vector = std::vector<const triangle*>;
    trigrid(const point3& min, const point3& max, const index3& dim);
    trigrid(const trigrid& grid);
    trigrid& operator=(const trigrid& rhs);
    ~trigrid();
    const triangle* insert(const triangle& tri);
    point3 center() const;
    point3 size() const;
    const index3& dim() const;
    point3 delta() const;
    void clear();
    bool empty() const;
    int count() const;
    int occupancy() const;
    const triangle_p_vector& operator[](const index3& i) const;
    void for_each_cell(std::function<void(const index3&, const triangle_p_vector&)> callback) const;
private:
    point3 _grid_min;
    point3 _grid_max;
    index3 _grid_dim;
    std::vector<triangle_p_vector> _grid_map;
    triangle_p_vector _triangle_p_vec;
};

#endif