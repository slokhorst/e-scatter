/**
 * @file src/cpl/xml.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "xml.h"
#include <fstream>
#include <sstream>
#include <cpl/text.h>

namespace cpl {
namespace xml {

element::element(const std::string& filename) {
	_name_tag = "xml::file";
	_attr_map["source"] = filename;
	std::stringstream ss;
	ss << text::file_to_string(filename);
	parser is(ss);
	while(parse_from(is));
}

element::element(const element& obj) {
	*this = obj;
}

element::~element() {
	clear();
}

element& element::operator=(const element& obj) {
	if(&obj == this)
		return *this;
	clear();
	_name_tag = obj._name_tag;
	_attr_map = obj._attr_map;
	for(element* child_p : obj._child_p_vec)
		_child_p_vec.push_back(new element(*child_p));
	return *this;
}

std::vector<const element*> element::children(const std::string& name_tag) const {
	std::vector<const element*> child_p_vec;
	for(const element* child_p : _child_p_vec)
		if(child_p->name_tag() == name_tag)
			child_p_vec.push_back(child_p);
	return child_p_vec;
}

void element::clear() {
	_name_tag.clear();
	_attr_map.clear();
	for(const element* child_p : _child_p_vec)
		delete child_p;
	_child_p_vec.clear();
}

bool element::parse_from(parser& is) {
	const std::string name_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789@-_";
	{/* parse text */
		const int row = is.row();
		const std::string xml_text = text::strip_string(is.get_while_not("<"), " \n\t");
		is.get();
		if(!xml_text.empty()) {
			xml::element* child_p = new xml::element;
			child_p->_name_tag = "xml::text";
			child_p->_attr_map["string"] = xml_text;
			child_p->_row = row;
			_child_p_vec.push_back(child_p);
		}
	}
	{/* parse child element */
		const int row = is.row();
		is.get_white_space();
		if(is.eof()) {
			if(_name_tag != "xml::file")
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
					"unexpected end of stream"
				);
			return false;
		} else if(is.peek() == '!') {
			/* parse comment */
			is.get();
			for(int i = 0; i < 2; i++)
				if(is.peek() == '-')
					is.get();
				else if(is.eof()) {
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"unexpected end of stream"
					);
				} else {
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
						"unexpected character `", static_cast<char>(is.peek()), "`",
						" in row ", row
					));
				}
			std::string xml_comment;
			while(true) {
				xml_comment += text::strip_string(is.get_while_not(">"), " \t\n");
				if(is.eof())
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"unexpected end of stream"
					);
				is.get();
				if(xml_comment.size() < 2)
					continue;
				if(*(xml_comment.crbegin()) == '-')
					if(*(xml_comment.crbegin()+1) == '-')
						break;
			}
		} else if(is.peek() == '/') {
			/* close tag */
			is.get();
			is.get_white_space();
			const std::string name_tag = text::strip_string(is.get_while(name_chars), " \t\n");
			if(is.peek() == '>')
				is.get();
			else if(is.eof()) {
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
					"unexpected end of stream"
				);
			} else {
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
					"unexpected character `", static_cast<char>(is.peek()), "`",
					" in row ", row
				));
			}
			if(_name_tag != name_tag) {
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
					"unexpected close tag `", name_tag, "`",
					" in row ", row,
					" while parsing element `", _name_tag, "`",
					" from row ", _row
				));
			}
			return false;
		} else if(isalpha(is.peek())) {
			/* open, empty or close tag */
			element* child_p = new element;
			child_p->_row = row;
			_child_p_vec.push_back(child_p);
			if(std::isalpha(is.peek()))
				child_p->_name_tag = text::strip_string(is.get_while(name_chars), " \t\n");
			else if(is.eof()) {
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
					"unexpected end of stream"
				);
			} else {
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
					"unexpected character `", static_cast<char>(is.peek()), "`",
					" in row ", row
				));
			}
			bool empty_tag = false;
			while(true) {
				is.get_white_space();
				if(is.eof()) {
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"unexpected end of stream"
					);
				} else if(is.peek() == '/') {
					if(empty_tag) {
						throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
							"duplicate symbol `", static_cast<char>(is.peek()), "`",
							" in row ", row
						));
					}
					is.get();
					empty_tag = true;
				} else if(is.peek() == '>') {
					is.get();
					if(!empty_tag)
						while(child_p->parse_from(is));
					break;
				} else if(isalnum(is.peek())) {
					const std::string attr_key = text::strip_string(is.get_while(name_chars), " \t\n");
					is.get_white_space();
					if(is.peek() == '=')
						is.get();
					else if(is.eof()) {
						throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
							"unexpected end of stream"
						);
					} else {
						throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
							"unexpected character `", static_cast<char>(is.peek()), "`"
							" in row ", row
						));
					}
					is.get_white_space();
					std::string symbol;
					if(std::string("\"'`").find_first_of(is.peek()) != std::string::npos)
						symbol.push_back(is.get());
					else if(is.eof()) {
						throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
							"unexpected end of stream"
						);
					} else {
						throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
							"unexpected character `", static_cast<char>(is.peek()), "`",
							" in row ", row
						));
					}
					const std::string attr_value = text::strip_string(is.get_while_not(symbol), " \t\n");
					if(is.eof())
						throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
							"unexpected end of stream"
						);
					is.get();
					child_p->_attr_map[attr_key] = attr_value;
				} else {
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
						"unexpected character `", static_cast<char>(is.peek()), "`",
						" in row ", row
					));
				}
			}
		} else {
			throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
				"unexpected character `", static_cast<char>(is.peek()), "`",
				" in row ", row
			));
		}
	}
	return true;
}

schema& schema::set_required_attribute(const std::string& key) {
	_required_attr_set.insert(key);
	return *this;
}

schema& schema::set_optional_attribute(const std::string& key) {
	_optional_attr_set.insert(key);
	return *this;
}

schema& schema::set_optional_child(const std::string& name_tag) {
	_optional_child_set.insert(name_tag);
	return *this;
}

void schema::apply_to(const element& parent) const {
	if(parent.name_tag() != "xml::file") {
		for(const std::string& attr_key : _required_attr_set)
			if(parent.attributes().count(attr_key) == 0)
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
					"missing required attribute `", attr_key, "`",
					" in element `", parent.name_tag(), "`",
					" in row ", parent.row()
				));
		for(auto cit = parent.attributes().cbegin(); cit != parent.attributes().cend(); cit++)
			if(_required_attr_set.count(cit->first) > 0)
				continue;
			else if(_optional_attr_set.count(cit->first) > 0)
				continue;
			else {
//				static_function_log.warning(text::cat(
//					"ignoring attribute `", cit->first, "`" 
//					" from element `", parent.name_tag(), "` in row ", parent.row()
//				));
			}
	}
	for(const element* child_p : parent.children())
		if(child_p->name_tag() == "xml::text") {
//			static_function_log.warning(text::cat(
//				"ignoring text from row ", child_p->row()
//			));
		} else if(_optional_child_set.count(child_p->name_tag()) == 0) {
//			static_function_log.warning(text::cat(
//				"ignoring element `", child_p->name_tag(), "` in row ", child_p->row()
//			));
		}
}

}
}