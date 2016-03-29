/**
 * @file src/cdsem/point3.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__POINT3__HEADER_INCLUDED
#define eSCATTER__CDSEM__POINT3__HEADER_INCLUDED

struct point3 {
    point3(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    point3() : point3(0, 0, 0) {};
    double x, y, z;
};

#endif