/**
 * @file src/edge-detect/matrix.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "matrix.h"

matrix::matrix(const size_t& num_rows, const size_t& num_cols) {
	_num_rows = num_rows;
	_num_cols = num_cols;
	_element_vec.assign(_num_rows*_num_cols, 0);
}

matrix& matrix::copy_from(const matrix& obj) {
	return *this = obj;
}

matrix matrix::transpose() const {
	matrix transposed_matrix(_num_cols, _num_rows);
	for(size_t j = 0; j < _num_cols; j++)
		for(size_t i = 0; i < _num_rows; i++)
			transposed_matrix.at(j,i) = at(i,j);
	return transposed_matrix;
}

const double& matrix::at(const size_t& i, const size_t& j) const {
	return _element_vec.at(i+j*_num_rows);
}

double& matrix::at(const size_t& i, const size_t& j) {
	return _element_vec.at(i+j*_num_rows);
}

const double* matrix::data() const {
	return _element_vec.data();
}

double* matrix::mutable_data() {
	return _element_vec.data();
}

const size_t& matrix::num_rows() const {
	return _num_rows;
}

const size_t& matrix::num_cols() const {
	return _num_cols;
}
