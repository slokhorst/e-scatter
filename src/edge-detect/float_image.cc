/**
 * @file src/edge-detect/float_image.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */


#include "float_image.h"
#include <Magick++.h>

float_image::float_image(const std::string& file) : m_image_matrix(1, 1) {
	const Magick::Image I(file);
	m_image_matrix = matrix(I.rows(), I.columns());
	for(size_t x = 0; x < width(); x++)
		for(size_t y = 0; y < height(); y++)
			(*this)(x,y) = Magick::ColorGray(I.pixelColor(x, y)).shade();
}

float_image float_image::transpose() const {
	float_image I(*this);
	I.m_image_matrix = I.m_image_matrix.transpose();
	return I;
}