/*!
 * @file src/cpl/random.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__RANDOM__HEADER_INCLUDED
#define CPL__RANDOM__HEADER_INCLUDED

#include <cinttypes>
#include <random>

namespace cpl {

/*!
 *
 */
class random {
public:
	inline random();
	random(const random&) = delete;
	~random() = default;
	random& operator=(const random&) = delete;
	inline double uniform();
	inline uint32_t poisson(double mu);
	inline double gaussian(double mu, double sigma);
private:
	std::random_device _random_device;
	std::mt19937_64 _random_generator;
	std::uniform_real_distribution<double> _uniform_distribution;
};

}

#include <cpl/random.inl>

#endif