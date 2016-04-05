/**
 * @file src/cpl/text.inl
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__TEXT__INCLUDE_INCLUDED
#define CPL__TEXT__INCLUDE_INCLUDED

#include <sstream>

namespace cpl {
namespace text {

inline std::string cat() {
	return std::string();
}

template<typename Arg, typename...Args> std::string cat(const Arg& x, const Args&... z) {
	std::ostringstream oss;
	oss << x << cat(z...);
	return oss.str();
}

}
}

#endif