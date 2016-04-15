#ifndef CPL__TRIANGLE__HEADER_INCLUDED
#define CPL__TRIANGLE__HEADER_INCLUDED

#include <array>
#include <cpl/vector3.h>

namespace cpl {

class triangle {
public:
	inline triangle() {
		const double h = std::sqrt(3.0/4);
		_vertex_array[0] = vector3(-0.5, -h/2, 0);
		_vertex_array[1] = vector3(0.5, -h/2, 0);
		_vertex_array[2] = vector3(0, h/2, 0);
	}
	inline triangle(const std::array<vector3,3>& vertices) {
		_vertex_array = vertices;
	}
	triangle(const triangle&) = default;
	~triangle() = default;
	triangle& operator=(const triangle&) = default;
	inline const std::array<vector3,3>& vertices() const {
		return _vertex_array;
	}
	inline std::array<vector3,3>& vertices() {
		return _vertex_array;
	}
	inline vector3 normal() const {
		const vector3 AB = _vertex_array[1]-_vertex_array[0];
		const vector3 AC = _vertex_array[2]-_vertex_array[0];
		return AB.cross_product(AC);
	}
private:
	std::array<vector3,3> _vertex_array;
};

}

#endif