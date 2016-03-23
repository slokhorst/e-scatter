/**
 * @file src/cdsem/trimesh.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "trimesh.hh"

bool trimesh::empty() const {
    return _triangle_vec.empty();
}

int trimesh::count() const {
    return _triangle_vec.size();
}

void trimesh::clear() {
    _triangle_vec.clear();
}

void trimesh::push(const triangle& tri) {
    _triangle_vec.push_back(tri);
}

void trimesh::push(const trimesh& mesh) {
    for(auto cit = mesh.cbegin(); cit != mesh.cend(); cit++)
        push(*cit);
}

const triangle& trimesh::operator[](const int i) const {
    return _triangle_vec.at(i);
}

triangle& trimesh::operator[](const int i) {
    return _triangle_vec.at(i);
}

std::vector<triangle>::iterator trimesh::begin() {
    return _triangle_vec.begin();
}

std::vector<triangle>::iterator trimesh::end() {
    return _triangle_vec.end();
}

std::vector<triangle>::const_iterator trimesh::cbegin() const {
    return _triangle_vec.cbegin();
}

std::vector<triangle>::const_iterator trimesh::cend() const {
    return _triangle_vec.cend();
}