/**
 * @file src/common/parser.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <compile-mat/parser.hh>
#include <common/constant.hh>

parser::parser() {
	_p.DefineConst("pi",constant::pi);
	_p.DefineConst("c",constant::c);
	_p.DefineConst("ec",constant::ec);
	_p.DefineConst("eV",constant::ec);
	_p.DefineConst("h",constant::h);
	_p.DefineConst("hbar",constant::hbar);
	_p.DefineConst("m",1);
	_p.DefineConst("sr",1);
}

double parser::eval(const std::string& expr) {
	_p.SetExpr(expr);
	return _p.Eval();
}
