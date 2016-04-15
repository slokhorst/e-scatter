#ifndef CPL__VECTOR3__HEADER_INCLUDED
#define CPL__VECTOR3__HEADER_INCLUDED

#include <cmath>

namespace cpl {

class vector3 {
public:
	vector3() = default;
	inline vector3(double _x, double _y, double _z) {
		x = _x;
		y = _y;
		z = _z;
	}
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
	inline vector3 rotate(const vector3& r) const;
	double x = 0;
	double y = 0;
	double z = 0;
};

inline vector3 operator*(double s, const vector3& r) {
	return r*s;
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

vector3 vector3::rotate(const vector3& t) const {
	const double theta = t.norm();
	if(theta == 0)
		return *this;
	const vector3 a = t/theta;
	const vector3 u = dot_product(a)*a;
	const vector3 v = (*this)-u;
	const vector3 w = a.cross_product(v);
	return u+v*std::cos(theta)+w*std::sin(theta);
}

}

#endif