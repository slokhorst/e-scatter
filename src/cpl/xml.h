/*!
 * @file src/cpl/xml.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__XML__HEADER_INCLUDED
#define CPL__XML__HEADER_INCLUDED

#include <map>
#include <set>
#include <string>
#include <vector>
#include <cpl/parser.h>

namespace cpl {
namespace xml {

class element {
public:
	element(const std::string& filename);
	element(const element&);
	~element();
	element& operator=(const element&);
	inline const std::string& name_tag() const;
	inline const std::map<std::string,std::string>& attributes() const;
	inline const std::vector<element*>& children() const;
	std::vector<const element*> children(const std::string&) const;
	inline const int& row() const;
private:
	element() = default;
	void clear();
	bool parse_from(parser&);
	std::string _name_tag;
	std::map<std::string,std::string> _attr_map;
	std::vector<element*> _child_p_vec;
	int _row = 1;
};

class schema {
public:
	schema() = default;
	schema(const schema&) = default;
	~schema() = default;
	schema& operator=(const schema&) = default;
	schema& copy_from(const schema&);
	schema& set_required_attribute(const std::string&);
	schema& set_optional_attribute(const std::string&);
	schema& set_optional_child(const std::string&);
	void apply_to(const element&) const;
private:
	std::set<std::string> _required_attr_set;
	std::set<std::string> _optional_attr_set;
	std::set<std::string> _optional_child_set;
};

}
}

#include <cpl/xml.inl>

#endif