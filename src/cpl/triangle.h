/*!
 * @file src/cpl/triangle.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__TRIANGLE__HEADER_INCLUDED
#define CPL__TRIANGLE__HEADER_INCLUDED

#include <array>
#include <cpl/vector3.h>

namespace cpl {

/*!
 *
 */
class triangle {
public:
	inline triangle();
	inline triangle(const std::array<vector3,3>& vertices);
	triangle(const triangle&) = default;
	~triangle() = default;
	triangle& operator=(const triangle&) = default;
	/*!
	 * @warning the result is always false.
	 */
	inline bool point_test(const vector3& p) const;
	inline vector3 centroid() const;
	vector3 min() const;
	vector3 max() const;
	triangle& rotate(const vector3& r);
	triangle& scale(const vector3& s);
	triangle& translate(const vector3& t);
	inline const std::array<vector3,3>& vertices() const;
	inline std::array<vector3,3>& vertices();
	/*!
	 * @warning the result is not normalized.
	 */
	inline vector3 normal() const;
private:
	std::array<vector3,3> _vertex_array;
};

}

#include <cpl/triangle.inl>

#endif