/**
 * @file src/cdsem/trimesh.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__TRIMESH__HEADER_INCLUDED
#define eSCATTER__CDSEM__TRIMESH__HEADER_INCLUDED

#include <vector>
#include "point3.hh"

class trimesh {
public:
    struct face {
        point3 A, B, C;
        int in, out;
    };
    bool empty() const;
    int size() const;
    void clear();
    void push(const face&);
    void push(const trimesh&);
    const face& operator[](const int) const;
    face& operator[](const int);
    std::vector<face>::iterator begin();
    std::vector<face>::iterator end();
    std::vector<face>::const_iterator cbegin() const;
    std::vector<face>::const_iterator cend() const;
private:
    std::vector<face> _face_vec;
};

#endif