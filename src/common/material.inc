/**
 * @file src/cdsem/material.inc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__MATERIAL__INLINE_INCLUDED
#define eSCATTER__CDSEM__MATERIAL__INLINE_INCLUDED

#include <cmath>
#include "interpolate.hh"

const std::string& material::name() const {
    return _name;
}

double material::fermi() const {
    return _fermi;
}

double material::barrier() const {
    return _barrier;
}

const optional<double>& material::band_gap() const {
    return _band_gap;
}

double material::phonon_loss() const {
    return _phonon_loss;
}

double material::density() const {
    return _density;
}

double material::elastic_tcs(double K) const {
    return std::exp(interpolate(_elastic_tcs, std::log(K)));
}

double material::elastic_icdf(double K, double P) const {
    return interpolate(_elastic_icdf, std::log(K), P);
}

double material::inelastic_tcs(double K) const {
    return std::exp(interpolate(_inelastic_tcs, std::log(K)));
}

double material::inelastic_icdf(double K, double P) const {
    return std::exp(interpolate(_inelastic_icdf, std::log(K), P));
}

#endif
