/**
 * @file src/cpl/vector3.inl
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__VECTOR3__INLINE_INCLUDED
#define CPL__VECTOR3__INLINE_INCLUDED

#include <cmath>

namespace cpl {

vector3::vector3(double _x, double _y, double _z) {
	x = _x;
	y = _y;
	z = _z;
}

vector3 vector3::operator-() const {
	return vector3(*this)*(-1);
}

vector3 vector3::operator*(double s) const {
	return vector3(s*x, s*y, s*z);
}

vector3 vector3::operator/(double s) const {
	return vector3(x/s, y/s, z/s);
}

vector3 vector3::operator+(const vector3& r) const {
	return vector3(x+r.x, y+r.y, z+r.z);
}

vector3 vector3::operator-(const vector3& r) const {
	return vector3(x-r.x, y-r.y, z-r.z);
}

vector3& vector3::operator*=(double s) {
	return *this = (*this)*s;
}

vector3& vector3::operator/=(double s) {
	return *this = (*this)/s;
}

vector3& vector3::operator+=(const vector3& r) {
	return *this = (*this)+r;
}

vector3& vector3::operator-=(const vector3& r) {
	return *this = (*this)-r;
}

double vector3::norm_pow_2() const {
	return dot_product(*this);
}

double vector3::norm() const {
	return std::sqrt(norm_pow_2());
}

double vector3::dot_product(const vector3& r) const {
	return x*r.x+y*r.y+z*r.z;
}

vector3 vector3::cross_product(const vector3& r) const {
	return vector3(r.z*y-r.y*z, r.x*z-r.z*x, r.y*x-r.x*y);
}

vector3 operator*(double s, const vector3& r) {
	return r*s;
}

}

#endif