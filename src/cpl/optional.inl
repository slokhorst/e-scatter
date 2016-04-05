/**
 * @file src/cpl/optional.inl
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__OPTIONAL__INLINE_INCLUDED
#define CPL__OPTIONAL__INLINE_INCLUDED

namespace cpl {

template<typename T>
optional<T>::optional(const T& value) {
	*this = value;
}

template<typename T>
optional<T>::optional(const optional& obj) {
	*this = obj;
}

template <typename T> optional<T>::~optional() {
	clear();
}

template<typename T>
optional<T>& optional<T>::operator=(const optional<T>& obj) {
	if(&obj == this)
		return *this;
	clear();
	if(obj.is_defined())
		*this = obj();
	return *this;
}

template<typename T>
optional<T>& optional<T>::operator=(const T& value) {
	if(_value_p == nullptr)
		_value_p = new T(value);
	else
		*_value_p = value;
	return *this;
}

template<typename T>
const T& optional<T>::operator()() const {
	if(!is_defined())
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ "value not defined");
	return *_value_p;
}

template<typename T>
bool optional<T>::is_defined() const {
	return (_value_p != nullptr);
}

template<typename T>
optional<T>& optional<T>::clear() {
	delete _value_p;
	_value_p = nullptr;
	return *this;
}

}

#endif