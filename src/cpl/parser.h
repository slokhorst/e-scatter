/*!
 * @file src/cpl/parser.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef CPL__PARSER__HEADER_INCLUDED
#define CPL__PARSER__HEADER_INCLUDED

#include <iostream>
#include <string>

namespace cpl {

/*!
 *
 */
class parser {
public:
	/*!
	 * @param[in,out] is
	 */
	parser(std::istream& is);
	parser(const parser&) = delete;
	~parser() = default;
	parser& operator=(const parser&) = delete;
	int row() const;
	bool eof();
	std::streambuf::int_type peek();
	std::streambuf::int_type get();
	/*!
	 * @param[in] chars
	 * @return
	 */
	std::string get_while(const std::string& chars);
	/*!
	 * @param[in] chars
	 * @return
	 */
	std::string get_while_not(const std::string& chars);
	/*!
	 * @return
	 */
	std::string get_white_space();
	/*!
	 * @return
	 */
	std::string get_digit_string();
	/*!
	 * @return
	 */
	std::string get_alpha_string();
	/*!
	 * @return
	 */
	std::string get_alnum_string();
	/*!
	 * @return
	 */
	std::string get_numeric();
private:
	std::streambuf* _streambuf_p;
	int _row = 1;
};

}

#endif