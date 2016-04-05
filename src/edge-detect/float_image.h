/**
 * @file src/edge-detect/float_image.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef EDGE_DETECT__FLOAT_IMAGE__HEADER_INCLUDED
#define EDGE_DETECT__FLOAT_IMAGE__HEADER_INCLUDED

#include <string>
#include "matrix.h"

class float_image {
public:
	float_image(const std::string& file);

	const double& operator()(const size_t& ix, const size_t& iy) const { return m_image_matrix.at(iy,ix); }
	double& operator()(const size_t& ix, const size_t& iy) { return m_image_matrix.at(iy,ix); }

	size_t width() const { return m_image_matrix.num_cols(); }
	size_t height() const { return m_image_matrix.num_rows(); }

	float_image transpose() const;

private:
	matrix m_image_matrix;
};

#endif