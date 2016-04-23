/**
 * @file src/common/parser.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__PARSER__HEADER_INCLUDED
#define eSCATTER__COMMON__PARSER__HEADER_INCLUDED

#include <muParser.h>

class parser {
public:
	parser();
	double eval(const std::string& expr);
private:
	mu::Parser _p;
};

#endif
