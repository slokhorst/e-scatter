/**
 * @file src/cpl/xml.inl
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__XML__INLINE_INCLUDED
#define CPL__XML__INLINE_INCLUDED

namespace cpl {
namespace xml {

const std::string& element::name_tag() const {
	return _name_tag;
}

const std::map<std::string,std::string>& element::attributes() const {
	return _attr_map;
}

const std::vector<element*>& element::children() const {
	return _child_p_vec;
}

const int& element::row() const {
	return _row;
}

}
}

#endif