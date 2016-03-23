/**
 * @file src/cdsem/trimesh.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__TRIMESH__HEADER_INCLUDED
#define eSCATTER__CDSEM__TRIMESH__HEADER_INCLUDED

#include <vector>
#include "point3.hh"
#include "triangle.hh"

class trimesh {
public:
    bool empty() const;
    int count() const;
    void clear();
    void push(const triangle&);
    void push(const trimesh&);
    const triangle& operator[](const int) const;
    triangle& operator[](const int);
    std::vector<triangle>::iterator begin();
    std::vector<triangle>::iterator end();
    std::vector<triangle>::const_iterator cbegin() const;
    std::vector<triangle>::const_iterator cend() const;
private:
    std::vector<triangle> _triangle_vec;
};

#endif