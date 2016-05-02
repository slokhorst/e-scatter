/**
 * @file src/cdsem/material.inc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__MATERIAL__INLINE_INCLUDED
#define eSCATTER__CDSEM__MATERIAL__INLINE_INCLUDED

#include <cmath>
#include <common/interpolate.hh>

const std::string& material::name() const {
    return _name;
}

double material::fermi() const {
    return _fermi;
}

double material::barrier() const {
    return _barrier;
}

const optional<double>& material::bandgap() const {
    return _bandgap;
}

double material::density() const {
    return _density;
}

double material::elastic_tcs(double K) const {
    return std::exp(interpolate(_elastic_tcs, std::log(K)));
}

double material::elastic_dcs(double K, double P) const {
    return interpolate(_elastic_dcs, std::log(K), P);
}

double material::inelastic_tcs(double K) const {
    return std::exp(interpolate(_inelastic_tcs, std::log(K)));
}

double material::inelastic_dcs(double K, double P) const {
    return interpolate(_inelastic_dcs, std::log(K), P);
}

#endif