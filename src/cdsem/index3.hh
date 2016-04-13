/**
 * @file src/cdsem/index3.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__INDEX3__HEADER_INCLUDED
#define eSCATTER__CDSEM__INDEX3__HEADER_INCLUDED

struct index3 {
    index3(int _x, int _y, int _z) : x(_x), y(_y), z(_z) {};
    index3() : index3(0, 0, 0) {};
    int x, y, z;
};

#endif