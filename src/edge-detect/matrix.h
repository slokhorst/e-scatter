/**
 * @file src/edge-detect/matrix.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CCPO__MATRIX__HEADER_INCLUDED
#define CCPO__MATRIX__HEADER_INCLUDED

#include <vector>
#include <cstdlib>

class matrix {
public:
	matrix(const size_t& num_rows, const size_t& num_cols);
	matrix(const matrix&) = default;
	matrix& operator=(const matrix&) = default;
	matrix& copy_from(const matrix&);
	const double& at(const size_t& i, const size_t& j) const;
	double& at(const size_t& i, const size_t& j);
	const double* data() const;
	double* mutable_data();
	const size_t& num_rows() const;
	const size_t& num_cols() const;
	matrix transpose() const;
private:
	size_t _num_rows = 0;
	size_t _num_cols = 0;
	std::vector<double> _element_vec;
};

#endif