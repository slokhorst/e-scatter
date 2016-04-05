/*!
 * @file src/cpl/triangle.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "triangle.h"

namespace cpl {

vector3 triangle::min() const {
	double x = std::min({_vertex_array[0].x, _vertex_array[1].x, _vertex_array[2].x});
	double y = std::min({_vertex_array[0].y, _vertex_array[1].y, _vertex_array[2].y});
	double z = std::min({_vertex_array[0].z, _vertex_array[1].z, _vertex_array[2].z});
	return vector3(x, y, z);
}

vector3 triangle::max() const {
	double x = std::max({_vertex_array[0].x, _vertex_array[1].x, _vertex_array[2].x});
	double y = std::max({_vertex_array[0].y, _vertex_array[1].y, _vertex_array[2].y});
	double z = std::max({_vertex_array[0].z, _vertex_array[1].z, _vertex_array[2].z});
	return vector3(x, y, z);
}

triangle& triangle::rotate(const vector3& r) {
	for(vector3& vertex : _vertex_array)
		vertex = vertex.rotate(r);
	return *this;
}

triangle& triangle::scale(const vector3& s) {
	for(vector3& vertex : _vertex_array) {
		vertex.x *= s.x;
		vertex.y *= s.y;
		vertex.z *= s.z;
	}
	return *this;
}

triangle& triangle::translate(const vector3& t) {
	for(vector3& vertex : _vertex_array) {
		vertex.x += t.x;
		vertex.y += t.y;
		vertex.z += t.z;
	}
	return *this;
}

}