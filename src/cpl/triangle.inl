/**
 * @file src/cpl/triangle.inl
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__TRIANGLE__INLINE_INCLUDED
#define CPL__TRIANGLE__INLINE_INCLUDED

#include <algorithm>
#include <cmath>

namespace cpl {

triangle::triangle() {
	const double h = std::sqrt(3.0/4);
	_vertex_array[0] = vector3(-0.5, -h/2, 0);
	_vertex_array[1] = vector3(0.5, -h/2, 0);
	_vertex_array[2] = vector3(0, h/2, 0);
}

triangle::triangle(const std::array<vector3,3>& vertices) {
	_vertex_array = vertices;
}

bool triangle::point_test(const vector3& p) const {
	return false;
}

vector3 triangle::centroid() const {
	return (_vertex_array[0]+_vertex_array[1]+_vertex_array[2])/3;
}

const std::array<vector3,3>& triangle::vertices() const {
	return _vertex_array;
}

std::array<vector3,3>& triangle::vertices() {
	return _vertex_array;
}

vector3 triangle::normal() const {
	const vector3 AB = _vertex_array[1]-_vertex_array[0];
	const vector3 AC = _vertex_array[2]-_vertex_array[0];
	return AB.cross_product(AC);
}

}

#endif