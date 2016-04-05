/**
 * @file src/edge-detect/gauss2_model.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "gauss2_model.h"
#include <cfloat>
#include <cmath>
#include <vector>
#include "nelder_mead.h"

gauss2_model::gauss2_model(const float_image& image) {
	std::vector<double> profile(image.width());
	for(size_t j = 0; j < image.height(); j++)
		for(size_t i = 0; i < image.width(); i++)
			profile[i] += image(i,j)/image.height();
	double max_amplitude = 0;
	for(size_t i = 0; i < profile.size(); i++)
		max_amplitude = std::max(max_amplitude, profile[i]);
	const double profile_size = static_cast<double>(profile.size()-1);
	const nelder_mead::params<6> simplex_min = {
		0, /* center */
		max_amplitude*0.9, /* amplitude */
		0, /* base.first */
		0, /* base.second */
		1, /* std.first */
		1 /* std.second */
	};
	const nelder_mead::params<6> simplex_max = {
		profile_size,	/* center */
		max_amplitude*1.1,/* amplitude */
		max_amplitude, /* base.first */
		max_amplitude, /* base.second */
		profile_size/2, /* std.first */
		profile_size/2 /* std.second */
	};
	auto chi_2_fun = [&profile](const nelder_mead::params<6>& params)->double {
		gauss2_model model;
		model.center = params[0];
		model.amplitude = params[1];
		model.base.first = params[2];
		model.base.second = params[3];
		model.std.first = params[4];
		model.std.second = params[5];
		double chi_2 = 0;
		for(size_t i = 0; i < profile.size(); i++)
			chi_2 += (profile[i]-model(i))*(profile[i]-model(i));
		return chi_2;
	};
	auto simplex = nelder_mead::randomize<6>(simplex_min, simplex_max);
	while(nelder_mead::iterate<6>(chi_2_fun, simplex) > 10*DBL_EPSILON)
		nelder_mead::bound<6>(simplex, simplex_min, simplex_max);
	chi_2 = chi_2_fun(simplex[0]);
	center = simplex[0][0];
	amplitude = simplex[0][1];
	base.first = simplex[0][2];
	base.second = simplex[0][3];
	std.first = simplex[0][4];
	std.second = simplex[0][5];
};

double gauss2_model::operator()(const double& x) const {
	double b, s;
	if(x < center) {
		b = base.first;
		s = std.first;
	} else {
		b = base.second;
		s = std.second;
	}
	return b + (amplitude-b)*exp(-0.5*(x-center)*(x-center)/(s*s));
}