/**
 * @file src/cdsem/trimesh.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "trimesh.h"

bool trimesh::empty() const {
	return _face_vec.empty();
}

int trimesh::size() const {
	return _face_vec.size();
}

void trimesh::clear() {
	_face_vec.clear();
}

void trimesh::push(const face& triangle) {
	_face_vec.push_back(triangle);
}

void trimesh::push(const trimesh& triangle_mesh) {
	for(auto cit = triangle_mesh.cbegin(); cit != triangle_mesh.cend(); cit++)
		push(*cit);
}

const trimesh::face& trimesh::operator[](const int i) const {
	return _face_vec.at(i);
}

trimesh::face& trimesh::operator[](const int i) {
	return _face_vec.at(i);
}

std::vector<trimesh::face>::iterator trimesh::begin() {
	return _face_vec.begin();
}

std::vector<trimesh::face>::iterator trimesh::end() {
	return _face_vec.end();
}

std::vector<trimesh::face>::const_iterator trimesh::cbegin() const {
	return _face_vec.cbegin();
}

std::vector<trimesh::face>::const_iterator trimesh::cend() const {
	return _face_vec.cend();
}