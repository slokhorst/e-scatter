/**
 * @file src/edge-detect/nelder_mead.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CCPO__NELDER_MEAD__HEADER_INCLUDED
#define CCPO__NELDER_MEAD__HEADER_INCLUDED

#include <array>
#include <cmath>
#include <functional>

namespace nelder_mead {
template<int N> using params = std::array<double,N>;
template<int N> using function = std::function<double(const nelder_mead::params<N>&)>;
template<int N> using simplex = std::array<params<N>,N+1>;
template<int N> simplex<N> randomize(const params<N>& min, const params<N>& max);
template<int N> void bound(simplex<N>& simplex, const params<N>& min, const params<N>& max);
template<int N> double iterate(const function<N>& fun, simplex<N>& simplex);
}

#include "nelder_mead.inl"

#endif