/**
 * @file src/cdsem/point3.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__POINT3__HEADER_INCLUDED
#define eSCATTER__CDSEM__POINT3__HEADER_INCLUDED

struct point3 {
    point3();
    point3(double _x, double _y, double _z);
    point3 operator+(const point3& p) const;
    point3 operator-(const point3& p) const;
    point3& operator+=(const point3& p);
    point3& operator-=(const point3& p);
    point3 operator*(double s) const;
    point3 operator/(double s) const;
    point3& operator*=(double s);
    point3& operator/=(double s);
    double norm() const;
    double x, y, z;
};

double dot_product(const point3& p1, const point3& p2);
point3 cross_product(const point3& p1, const point3& p2);

#endif