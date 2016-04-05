/*!
 * @file src/cpl/vector3.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__VECTOR3__HEADER_INCLUDED
#define CPL__VECTOR3__HEADER_INCLUDED

namespace cpl {

/*!
 * A 3-dimensional Euclidean vector.
 */
class vector3 {
public:
	vector3() = default;
	inline vector3(double x, double y, double z);
	vector3(const vector3&) = default;
	~vector3() = default;
	vector3& operator=(const vector3&) = default;
	inline vector3 operator-() const;
	inline vector3 operator*(double x) const;
	inline vector3 operator/(double x) const;
	inline vector3 operator+(const vector3& r) const;
	inline vector3 operator-(const vector3& r) const;
	inline vector3& operator*=(double x);
	inline vector3& operator/=(double x);
	inline vector3& operator+=(const vector3& r);
	inline vector3& operator-=(const vector3& r);
	inline double norm_pow_2() const;
	inline double norm() const;
	inline double dot_product(const vector3& r) const;
	inline vector3 cross_product(const vector3& r) const;
	vector3 rotate(const vector3& r) const;
	double x = 0;
	double y = 0;
	double z = 0;
};

inline vector3 operator*(double x, const vector3& r);

}

#include <cpl/vector3.inl>

#endif