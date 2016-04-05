/**
 * @file src/cpl/optional.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__OPTIONAL__HEADER_INCLUDED
#define CPL__OPTIONAL__HEADER_INCLUDED

namespace cpl {

/*!
 *
 */
template<typename T>
class optional {
public:
	optional() = default;
	optional(const T&);
	optional(const optional&);
	~optional();
	optional& operator=(const optional&);
	optional& operator=(const T&);
	/*!
 	 * @warning throws an exception if the value is not defined.
 	 */
	const T& operator()() const;
	/*!
 	 * @return
 	 */
	bool is_defined() const;
	/*!
 	 * @return
 	 */
	optional& clear();
private:
	T* _value_p = nullptr;
};

}

#include <cpl/optional.inl>

#endif