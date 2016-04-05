/**
 * @file src/cpl/parser.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "parser.h"
#include <cpl/text.h>

namespace cpl {

parser::parser(std::istream& is) {
	_streambuf_p = is.rdbuf();
}

int parser::row() const {
	return _row;
}

bool parser::eof() {
	if(_streambuf_p->sgetc() == std::streambuf::traits_type::eof())
		return true;
	return false;
}

std::streambuf::int_type parser::peek() {
	return _streambuf_p->sgetc();
}

std::streambuf::int_type parser::get() {
	const std::streambuf::int_type c = _streambuf_p->sbumpc();
	if(std::streambuf::traits_type::eq_int_type(c, std::streambuf::traits_type::to_int_type('\n')))
		_row++;
	return c;
}

std::string parser::get_while(const std::string& c_str) {
	std::string str;
	while(c_str.find(peek()) != std::string::npos)
		str += get();
	return str;
}

std::string parser::get_while_not(const std::string& c_str) {
	std::string str;
	while(!eof() && (c_str.find(peek()) == std::string::npos))
		str += get();
	return str;
}

std::string parser::get_white_space() {
	std::string str;
	while(std::isspace(peek()))
		str += get();
	return str;
}

std::string parser::get_digit_string() {
	std::string str;
	while(std::isdigit(peek()))
		str += get();
	return str;
}

std::string parser::get_alpha_string() {
	std::string str;
	while(std::isalpha(peek()))
		str += get();
	return str;
}

std::string parser::get_alnum_string() {
	std::string str;
	while(std::isalnum(peek()))
		str += get();
	return str;
}

std::string parser::get_numeric() {
	std::string str;
	get_white_space();
	if((peek() == '+') || (peek() == '-'))
		str += get();
	str += get_digit_string();
	if(peek() == '.') {
		str += get();
		if(!std::isdigit(peek())) {
			if(eof())
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
					"unexpected end of stream"
				);
			else
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
					"unexpected character `", static_cast<char>(peek()), "`"
				));
		}
		str += get_digit_string();
	}
	if(str.empty())
		return str;
	if((peek() == 'e') || (peek() == 'E'))
		str += get();
	else {
		if((str == "+") || (str == "-"))
			throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
				"failed to interpret `", str, "` numerically"
			));
		return str;
	}
	if((peek() == '+') || (peek() == '-'))
		str += get();
	if(!std::isdigit(peek())) {
		if(eof())
			throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
				"unexpected end of stream"
			);
		else
			throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
				"unexpected character `", static_cast<char>(peek()), "`"
			));
	}
	str += get_digit_string();
	return str;
}

}