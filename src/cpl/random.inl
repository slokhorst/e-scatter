/**
 * @file src/cpl/random.inl
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__RANDOM__INLINE_INCLUDED
#define CPL__RANDOM__INLINE_INCLUDED

namespace cpl {

random::random() : _random_generator(_random_device()), _uniform_distribution(0, 1) {
}

double random::uniform() {
	return _uniform_distribution(_random_generator);
}

uint32_t random::poisson(double mu) {
	return std::poisson_distribution<uint32_t>(mu)(_random_generator);
}

double random::gaussian(double mu, double sigma) {
	return std::normal_distribution<double>(mu, sigma)(_random_generator);
}

}

#endif