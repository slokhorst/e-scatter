/**
 * @file src/cdsem/triangle.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__TRIANGLE__HEADER_INCLUDED
#define eSCATTER__CDSEM__TRIANGLE__HEADER_INCLUDED

#include "point3.hh"

struct triangle {
    enum id_enum : int {
        NOP = -128,
        TERMINATOR = -127,
        DETECTOR = -126,
        DETECTOR_LT50 = -125,
        DETECTOR_GE50 = -124,
        VACUUM = -123,
        MIRROR = -122
    };
    triangle(const point3& _A, const point3& _B, const point3& _C) : A(_A), B(_B), C(_C), in(NOP), out(NOP) {};
    triangle() : triangle(point3(0, 0, 0), point3(1, 0, 0), point3(0, 1, 0)) {};
    point3 A, B, C;
    int in, out;
};

#endif