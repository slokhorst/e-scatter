/**
 * @file src/edge-detect/gauss2_model.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef EDGE_DETECT__GAUSS2_MODEL__HEADER_INCLUDED
#define EDGE_DETECT__GAUSS2_MODEL__HEADER_INCLUDED

#include "float_image.h"

struct gauss2_model {
	gauss2_model() = default;
	gauss2_model(const float_image& image);
	double operator()(const double& x) const;

	double chi_2 = 0;
	double center = 0;
	double amplitude = 1;
	std::pair<double, double> base = { 0, 0 };
	std::pair<double, double> std = { 1, 1 };
};

#endif