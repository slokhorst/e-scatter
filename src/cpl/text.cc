/**
 * @file src/cpl/text.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "text.h"
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <cpl/parser.h>

namespace cpl {
namespace text {

std::string file_to_string(const std::string& file) {
	std::string str;
	std::ifstream ifs(file);
	if(!ifs.good())
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
			"failed to open file `", file, "` for reading"
		));
	while(true) {
		const auto c = ifs.get();
		if(ifs.eof())
			break;
		else if(ifs.bad() || ifs.fail())
			throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
				"failed to read from file `", file, "`"
			));
		str.push_back(c);
	}
	ifs.close();
	return str;
}

std::string strip_string(const std::string& str, const std::string& c_str) {
	const auto i = str.find_first_not_of(c_str);
	if(i == std::string::npos)
		return std::string();
	const auto j = str.find_last_not_of(c_str);
	return str.substr(i,j-i+1);
}

std::vector<std::string> split_string(const std::string& str, char c) {
	std::vector<std::string> str_vec;
	std::string sub_str;
	for(auto cit = str.cbegin(); cit < str.cend(); cit++)
		if(*cit != c)
			sub_str.push_back(*cit);
		else {
			str_vec.push_back(sub_str);
			sub_str.clear();
		}
	if(!sub_str.empty())
		str_vec.push_back(sub_str);
	return str_vec;
}

std::string merge_string(const std::vector<std::string>& str_vec, char c) {
	if(str_vec.empty())
		return std::string();
	std::ostringstream oss;
	auto cit = str_vec.cbegin();
	oss << *cit;
	while(++cit != str_vec.cend())
		oss << c << *cit;
	return oss.str();
}

std::string space_to_tab(const std::string& str, uint8_t tab_size) {
	std::string lhs;
	int space_count = 0;
	for(auto cit = str.cbegin(); cit != str.cend(); cit++)
		switch(*cit) {
			case ' ': {
				if(++space_count == tab_size) {
					lhs.push_back('\t');
					space_count = 0;
				}
				break;
			}
			case '\t': {
				lhs.push_back('\t');
				space_count = 0;
				break;
			}
			default: {
				for(int i = 0; i < space_count; i++)
					lhs.push_back(' ');
				space_count = 0;
				lhs.push_back(*cit);
				break;
			}
		}
	for(int i = 0; i < space_count; i++)
		lhs.push_back(' ');
	return lhs;
}

std::string tab_to_space(const std::string& str, uint8_t tab_size) {
	std::string lhs;
	int space_count = 0;
	for(auto cit = str.cbegin(); cit != str.cend(); cit++)
		switch(*cit) {
			case ' ': {
				if(++space_count == tab_size)
					space_count = 0;
				lhs.push_back(' ');
				break;
			}
			case '\t': {
				for(int i = space_count; i < tab_size; i++)
					lhs.push_back(' ');
				space_count = 0;
				break;
			}
			default: {
				lhs.push_back(*cit);
				space_count = 0;
				break;
			}
		}
	return lhs;
}

std::string float64(double x) {
	std::ostringstream oss;
	oss << std::setprecision(DBL_DIG) << std::scientific << std::showpos << x;
	return oss.str();
}

double float64(const std::string& str) {
	std::stringstream ss(strip_string(str));
	double x;
	try {
		x = stod(parser(ss).get_numeric());
		ss.get();
		if(!ss.eof())
			throw std::exception();
	} catch(...) {
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
			"failed to interpret `", str, "`"
		));
	}
	return x;
}

std::string int32(int32_t x) {
	std::ostringstream oss;
	oss << x;
	return oss.str();
}

int32_t int32(const std::string& str) {
	return static_cast<int32_t>(float64(str));
}

std::string strc(const std::string& str, color_enum c) {
	std::map<color_enum, std::string> c_map = {
		{DEFAULT, "\033[0m"},
		{BLACK, "\033[;30m"},
		{RED, "\033[;31m"},
		{GREEN, "\033[;32m"},
		{YELLOW, "\033[;33m"},
		{BLUE, "\033[;34m"},
		{MAGENTA, "\033[;35m"},
		{CYAN, "\033[;36m"},
		{LIGHTGRAY, "\033[;37m"}
	};
	return c_map[c]+str+c_map[DEFAULT];
}

std::string strb(const std::string& str, color_enum c) {
	return strc(std::string("\033[1m")+str, c);
}

}
}