/*!
 * @file src/cpl/mathexpr.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__MATHEXPR__HEADER_INCLUDED
#define CPL__MATHEXPR__HEADER_INCLUDED

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <cpl/text.h>

namespace cpl {

/*!
 * A mathematical expression, which can contain variables
 * and can be evaluated to a number.
 * Has an accompanying textual human-readable representation,
 * which can be read by this class.
 */
class mathexpr {
public:
	mathexpr(std::istream& is) : _postfix_p_list(tokenize_stream(is)) {}
	mathexpr(const std::string& s);
	mathexpr(double x) : mathexpr(text::float64(x)) {}
	mathexpr(const mathexpr& obj) { *this = obj; }
	~mathexpr() { delete_tokens(); }
	mathexpr& operator=(const mathexpr& obj);
	std::string debug_string() const;
	double evaluate() const;
	mathexpr& substitute(const std::string& var, const mathexpr& expr);
	mathexpr& substitute(const std::string& var, double x);
	mathexpr& substitute(const std::map<std::string, mathexpr>& subs_map);
	mathexpr& substitute(const std::map<std::string, double>& subs_map);
private:
	class token;
	class numeric_token;
	class variable_token;
	class operator_token;
	class symbol_token;
	class function_token;
	std::list<token*> tokenize_stream(std::istream& is) const;
	std::list<token*> infix_to_postfix(std::list<token*>& infix_p_list) const;
	void delete_tokens();
	std::list<token*> _postfix_p_list;
};

}

#endif