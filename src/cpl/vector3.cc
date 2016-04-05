/**
 * @file src/cpl/vector3.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "vector3.h"

namespace cpl {

vector3 vector3::rotate(const vector3& t) const {
	const double theta = t.norm();
	if(theta == 0)
		return *this;
	const vector3 a = t/theta;
	const vector3 u = dot_product(a)*a;
	const vector3 v = (*this)-u;
	const vector3 w = a.cross_product(v);
	return u+v*std::cos(theta)+w*std::sin(theta);
}

}