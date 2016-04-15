/**
 * @file src/cdsem/point3.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#include "point3.hh"
#include <cmath>

point3::point3() : point3(0, 0, 0) {
};

point3::point3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {
}

point3 point3::operator+(const point3& p) const {
    return point3(x+p.x, y+p.y, z+p.z);
}

point3 point3::operator-(const point3& p) const {
    return point3(x-p.x, y-p.y, z-p.z);
}

point3& point3::operator+=(const point3& p) {
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
}

point3& point3::operator-=(const point3& p) {
    x -= p.x;
    y -= p.y;
    z -= p.z;
    return *this;
}

point3 point3::operator*(double s) const {
    return point3(x*s, y*s, z*s);
}

point3 point3::operator/(double s) const {
    return point3(x/s, y/s, z/s);
}

point3& point3::operator*=(double s) {
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

point3& point3::operator/=(double s) {
    x /= s;
    y /= s;
    z /= s;
    return *this;
}

double point3::norm() const {
    return std::sqrt(x*x+y*y+z*z);
}

double dot_product(const point3& p1, const point3& p2) {
    return p1.x*p2.x+p1.y*p2.y+p1.z*p2.z;
}

point3 cross_product(const point3& p1, const point3& p2) {
    return point3(p1.y*p2.z-p1.z*p2.y, p1.z*p2.x-p1.x*p2.z, p1.x*p2.y-p1.y*p2.x);
}