/**
 * @file src/cpl/mathexpr.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */
 
#include "mathexpr.h"
#include <cmath>
#include <limits>
#include <sstream>
#include <stack>
#include <common/constant.hh>
#include <cpl/parser.h>

namespace cpl {

enum token_enum {
	NUMERIC_TOKEN, VARIABLE_TOKEN, OPERATOR_TOKEN, SYMBOL_TOKEN, FUNCTION_TOKEN
};
enum operator_enum {
	ADD_OPERATOR, SUB_OPERATOR, MUL_OPERATOR, DIV_OPERATOR, NEG_OPERATOR, POW_OPERATOR
};
enum symbol_enum {
	OPEN_SYMBOL, CLOSE_SYMBOL
};
enum function_enum {
	SQRT_FUNCTION, EXP_FUNCTION, LOG_FUNCTION, SIN_FUNCTION, COS_FUNCTION, TAN_FUNCTION
};

class mathexpr::token {
public:
	virtual ~token() = default;
	virtual token_enum token_id() const = 0;
	virtual int precedence() const = 0;
	virtual std::string debug_string() const = 0;
	static token* copy(const token* token_p);
	static token* parse_stream(std::istream& is);
	static bool check(const token* token_p, const token_enum& id);
	static bool check(const token* token_p, const operator_enum& id);
	static bool check(const token* token_p, const symbol_enum& id);
	static bool check(const token* token_p, const function_enum& id);
};

class mathexpr::numeric_token final : public mathexpr::token {
public:
	numeric_token(double x) : m_float_value(x) {}
	token_enum token_id() const override final { return NUMERIC_TOKEN; }
	int precedence() const override final { return 0; }
	std::string debug_string() const override final { return text::cat("num:", text::float64(m_float_value)); }
	double value() const { return m_float_value; }

private:
	double m_float_value;
};

class mathexpr::variable_token final : public mathexpr::token {
public:
	variable_token(const std::string& name) : m_var_name(name) {}
	token_enum token_id() const override final { return VARIABLE_TOKEN; }
	int precedence() const override final { return 0; }
	std::string debug_string() const override final { return text::cat("var:", m_var_name); }
	bool is_equal_to(const token* token_p) const {
		if(token::check(token_p, VARIABLE_TOKEN)) {
			const auto* variable_p = dynamic_cast<const variable_token*>(token_p);
			if(m_var_name == variable_p->m_var_name)
				return true;
		}
		return false;
	}
	const std::string& name() const { return m_var_name; }

private:
	std::string m_var_name;
};

class mathexpr::operator_token final : public mathexpr::token {
public:
	operator_token(const operator_enum& id) : m_operator_id(id) {}
	token_enum token_id() const override final { return OPERATOR_TOKEN; }
	int precedence() const override final {
		switch(m_operator_id) {
			case ADD_OPERATOR: return -1;
			case SUB_OPERATOR: return -1;
			case MUL_OPERATOR: return -2;
			case DIV_OPERATOR: return -2;
			case NEG_OPERATOR: return +3;
			case POW_OPERATOR: return +3;
		}
		return 0;
	};
	std::string debug_string() const override final {
		switch(m_operator_id) {
			case ADD_OPERATOR: return std::string("op:+");
			case SUB_OPERATOR: return std::string("op:-");
			case MUL_OPERATOR: return std::string("op:*");
			case DIV_OPERATOR: return std::string("op:/");
			case NEG_OPERATOR: return std::string("op:~");
			case POW_OPERATOR: return std::string("op:^");
		}
		return std::string("op:?");
	}
	operator_enum operator_id() const { return m_operator_id; }

private:
	operator_enum m_operator_id;
};

class mathexpr::symbol_token final : public mathexpr::token {
public:
	symbol_token(const symbol_enum& id) : m_symbol_id(id) {}
	token_enum token_id() const override final { return SYMBOL_TOKEN; }
	int precedence() const override final { return 0; }
	std::string debug_string() const override final {
		switch(m_symbol_id) {
			case OPEN_SYMBOL: return std::string("sym:(");
			case CLOSE_SYMBOL: return std::string("sym:)");
		}
		return std::string("sym:?");
	}
	symbol_enum symbol_id() const { return m_symbol_id; }

private:
	symbol_enum m_symbol_id;
};

class mathexpr::function_token final : public mathexpr::token {
public:
	function_token(const function_enum& id) : m_function_id(id) {}
	token_enum token_id() const override final { return FUNCTION_TOKEN; }
	int precedence() const override final { return std::numeric_limits<int>::max(); }
	std::string debug_string() const override final {
		switch(m_function_id) {
			case SQRT_FUNCTION: return std::string("fun:sqrt");
			case EXP_FUNCTION: return std::string("fun:exp");
			case LOG_FUNCTION: return std::string("fun:log");
			case SIN_FUNCTION: return std::string("fun:sin");
			case COS_FUNCTION: return std::string("fun:cos");
			case TAN_FUNCTION: return std::string("fun:tan");
		}
		return std::string("fun:?");
	}
	function_enum function_id() const { return m_function_id; }

private:
	function_enum m_function_id;
};

mathexpr::token* mathexpr::token::copy(const token* token_p) {
	switch(token_p->token_id()) {
		case NUMERIC_TOKEN: {
			const auto* numeric_p = dynamic_cast<const numeric_token*>(token_p);
			return new numeric_token(*numeric_p);
		}
		case VARIABLE_TOKEN: {
			const auto* variable_p = dynamic_cast<const variable_token*>(token_p);
			return new variable_token(*variable_p);
		}
		case OPERATOR_TOKEN: {
			const auto* operator_p = dynamic_cast<const operator_token*>(token_p);
			return new operator_token(*operator_p);
			break;
		}
		case SYMBOL_TOKEN: {
			const auto* symbol_p = dynamic_cast<const symbol_token*>(token_p);
			return new symbol_token(*symbol_p);
			break;
		}
		case FUNCTION_TOKEN: {
			const auto* function_p = dynamic_cast<const function_token*>(token_p);
			return new function_token(*function_p);
		}
	}
	return nullptr;
}

mathexpr::token* mathexpr::token::parse_stream(std::istream& is) {
	const std::map<char, operator_enum> operator_map = {
		{'+', ADD_OPERATOR},
		{'-', SUB_OPERATOR},
		{'*', MUL_OPERATOR},
		{'/', DIV_OPERATOR},
		{'~', NEG_OPERATOR},
		{'^', POW_OPERATOR}
	};
	const std::map<char, symbol_enum> symbol_map = {
		{'(', OPEN_SYMBOL},
		{')', CLOSE_SYMBOL}
	};
	const std::map<std::string, function_enum> function_map = {
		{"sqrt", SQRT_FUNCTION},
		{"exp",  EXP_FUNCTION},
		{"log",  LOG_FUNCTION},
		{"sin",  SIN_FUNCTION},
		{"cos",  COS_FUNCTION},
		{"tan",  TAN_FUNCTION}
	};
	token* token_p = nullptr;
	parser(is).get_white_space();
	const char c = is.peek();
	if(operator_map.count(c) != 0) {
		const operator_enum operator_id = operator_map.find(c)->second;
		token_p = new operator_token(operator_id);
		is.get();
	} else if(symbol_map.count(c) != 0) {
		const symbol_enum symbol_id = symbol_map.find(c)->second;
		token_p = new symbol_token(symbol_id);
		is.get();
	} else {
		const std::string float_str = parser(is).get_numeric();
		if(!float_str.empty()) {
			const double value = text::float64(float_str);
			token_p = new numeric_token(value);
		} else {
			const std::string valid_name_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ@_:";
			const std::string name_str = parser(is).get_while(valid_name_chars);
			if(!name_str.empty()) {
				if(function_map.count(name_str) != 0) {
					const function_enum function_id = function_map.find(name_str)->second;
					token_p = new function_token(function_id);
				} else {
					token_p = new variable_token(name_str);
				}
			} else if(!is.eof()) {
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
					"unexpected character '", static_cast<char>(is.get()), "' detected while parsing token"
				));
			}
		}
	}
	return token_p;
}

bool mathexpr::token::check(const token* token_p, const token_enum& id) {
	if(token_p == nullptr)
		return false;
	if(token_p->token_id() == id)
		return true;
	return false;
}

bool mathexpr::token::check(const token* token_p, const operator_enum& id) {
	if(!token::check(token_p, OPERATOR_TOKEN))
		return false;
	if(dynamic_cast<const operator_token*>(token_p)->operator_id() == id)
		return true;
	return false;
}

bool mathexpr::token::check(const token* token_p, const symbol_enum& id) {
	if(!token::check(token_p, SYMBOL_TOKEN))
		return false;
	if(dynamic_cast<const symbol_token*>(token_p)->symbol_id() == id)
		return true;
	return false;
}

bool mathexpr::token::check(const token* token_p, const function_enum& id) {
	if(!token::check(token_p, FUNCTION_TOKEN))
		return false;
	if(dynamic_cast<const function_token*>(token_p)->function_id() == id)
		return true;
	return false;
}

mathexpr::mathexpr(const std::string& s) {
	std::stringstream ss(s);
	_postfix_p_list = tokenize_stream(ss);
}

mathexpr& mathexpr::operator=(const mathexpr& obj) {
	if(&obj == this)
		return *this;
	delete_tokens();
	for(const token* token_p : obj._postfix_p_list)
		_postfix_p_list.push_back(token::copy(token_p));
	return *this;
}

std::string mathexpr::debug_string() const {
	std::ostringstream oss;
	oss << "tokens: {";
	auto& token_p_list = _postfix_p_list;
	for(auto cit = token_p_list.cbegin(); cit != token_p_list.cend(); cit++) {
		if(cit == token_p_list.cbegin())
			oss << "'";
		else
			oss << ", '";
		oss << (*cit)->debug_string();
		oss << "'";
	}
	oss << "}";
	return oss.str();
}

double mathexpr::evaluate() const {
	const std::map<std::string, double> math_constant_map = {
		{"pi", constant::pi}
	};
	const std::map<std::string, double> physical_constant_map = {
		{"c", constant::c},
		{"ec", constant::ec},
		{"eps0", constant::eps0},
		{"eV", constant::ec},
		{"h", constant::h},
		{"hbar", constant::hbar},
		{"k", constant::k},
		{"me", constant::me},
		{"mu0", constant::mu0},
	};
	const std::map<std::string, double> scale_constant_map = {
		{ "m", 1}, {"mm", 1e-3}, {"um", 1e-6}, {"nm", 1e-9}, {"cm", 1e-2}, {"A", 1e-2},
		{ "s", 1}, {"ms", 1e-3}, {"us", 1e-6}, {"ns", 1e-9},
		{ "C", 1}, {"mC", 1e-3}, {"uC", 1e-6}, {"nC", 1e-9},
		{"kg", 1}, { "g", 1e-3}, {"mg", 1e-6}, {"ug", 1e-9},
		{"sr", 1}, {"mol", constant::NA}
	};
	const mathexpr expr(mathexpr(*this)
		.substitute(math_constant_map)
		.substitute(physical_constant_map)
		.substitute(scale_constant_map));
	std::stack<double> operand_stack;
	for(auto cit = expr._postfix_p_list.cbegin(); cit != expr._postfix_p_list.cend(); cit++)
		switch((*cit)->token_id()) {
			case NUMERIC_TOKEN: {
				const auto* numeric_p = dynamic_cast<const numeric_token*>(*cit);
				operand_stack.push(numeric_p->value());
				break;
			}
			case VARIABLE_TOKEN: {
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ text::cat(
					"undefined variable (", (*cit)->debug_string(), ") in postfix expression"
				));
			}
			case SYMBOL_TOKEN: {
				throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
					"symbol tokens are not allowed in postfix expressions"
				);
			}
			case OPERATOR_TOKEN: {
				const auto* operator_p = dynamic_cast<const operator_token*>(*cit);
				switch(operator_p->operator_id()) {
					case ADD_OPERATOR: {
						if(operand_stack.size() < 2)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"add operator requires two operands"
							);
						const double y = operand_stack.top();
						operand_stack.pop();
						const double x = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(x+y);
						break;
					}
					case SUB_OPERATOR: {
						if(operand_stack.size() < 2)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"sub operator requires two operands"
							);
						const double y = operand_stack.top();
						operand_stack.pop();
						const double x = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(x-y);
						break;
					}
					case MUL_OPERATOR: {
						if(operand_stack.size() < 2)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"mul operator requires two operands"
							);
						const double y = operand_stack.top();
						operand_stack.pop();
						const double x = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(x*y);
						break;
					}
					case DIV_OPERATOR: {
						if(operand_stack.size() < 2)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"div operator requires two operands"
							);
						const double y = operand_stack.top();
						operand_stack.pop();
						const double x = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(x/y);
						break;
					}
					case NEG_OPERATOR: {
						if(operand_stack.size() < 1)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"neg operator requires one operand"
							);
						const double z = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(-z);
						break;
					}
					case POW_OPERATOR: {
						if(operand_stack.size() < 2)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"pow operator requires two operands"
							);
						const double y = operand_stack.top();
						operand_stack.pop();
						const double x = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(std::pow(x,y));
						break;
					}
				}
				break;
			}
			case FUNCTION_TOKEN: {
				const auto* function_p = dynamic_cast<const function_token*>(*cit);
				switch(function_p->function_id()) {
					case SQRT_FUNCTION: {
						if(operand_stack.size() < 1)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"sqrt function requires one operand"
							);
						const double z = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(std::sqrt(z));
						break;
					}
					case EXP_FUNCTION: {
						if(operand_stack.size() < 1)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"exp function requires one operand"
							);
						const double z = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(std::exp(z));
						break;
					}
					case LOG_FUNCTION: {
						if(operand_stack.size() < 1)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"log function requires one operand"
							);
						const double z = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(std::log(z));
						break;
					}
					case SIN_FUNCTION: {
						if(operand_stack.size() < 1)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"sin function requires one operand"
							);
						const double z = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(std::sin(z));
						break;
					}
					case COS_FUNCTION: {
						if(operand_stack.size() < 1)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"cos function requires one operand"
							);
						const double z = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(std::cos(z));
						break;
					}
					case TAN_FUNCTION: {
						if(operand_stack.size() < 1)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"tan function requires one operand"
							);
						const double z = operand_stack.top();
						operand_stack.pop();
						operand_stack.push(std::tan(z));
						break;
					}
				}
				break;
			}
		}
	if(operand_stack.size() != 1)
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ "missing operators");
	if(!std::isfinite(operand_stack.top()))
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ "result is not finite");
	return operand_stack.top();
}

mathexpr& mathexpr::substitute(const std::string& var, const mathexpr& expr) {
	for(auto it = _postfix_p_list.begin(); it != _postfix_p_list.end(); it++) {
		if(!variable_token(var).is_equal_to(*it))
			continue;
		delete *it;
		it = _postfix_p_list.erase(it);
		std::vector<token*> token_p_vec;
		for(auto cit = expr._postfix_p_list.cbegin(); cit != expr._postfix_p_list.cend(); cit++)
			token_p_vec.push_back(token::copy(*cit));
		_postfix_p_list.insert(it, token_p_vec.cbegin(), token_p_vec.cend());
	}
	return *this;
}

mathexpr& mathexpr::substitute(const std::string& var, double x) {
	for(token*& token_p : _postfix_p_list)
		if(variable_token(var).is_equal_to(token_p)) {
			delete token_p;
			token_p = new numeric_token(x);
		}
	return *this;
}

mathexpr& mathexpr::substitute(const std::map<std::string, mathexpr>& subs_map) {
	for(const auto& map_item : subs_map)
		substitute(map_item.first, map_item.second);
	return *this;
}

mathexpr& mathexpr::substitute(const std::map<std::string, double>& subs_map) {
	for(const auto& map_item : subs_map)
		substitute(map_item.first, map_item.second);
	return *this;
}

std::list<mathexpr::token*> mathexpr::tokenize_stream(std::istream& is) const {
	std::list<token*> infix_p_list;
	int nesting_count = 0;
	while(true) {
		token* token_p = token::parse_stream(is);
		if(token_p == nullptr)
			break;
		switch(token_p->token_id()) {
			case NUMERIC_TOKEN: {
				if(infix_p_list.empty())
					break;
				if(token::check(infix_p_list.back(), NUMERIC_TOKEN))
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"consecutive numerical tokens detected"
					);
				if(token::check(infix_p_list.back(), VARIABLE_TOKEN))
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"numerical token after variable token detected"
					);
				break;
			}
			case VARIABLE_TOKEN: {
				if(infix_p_list.empty())
					break;
				if(token::check(infix_p_list.back(), NUMERIC_TOKEN))
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"variable token after numerical token detected"
					);
				if(token::check(infix_p_list.back(), VARIABLE_TOKEN))
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"consecutive variable tokens detected"
					);
				break;
			}
			case OPERATOR_TOKEN: {
				switch(dynamic_cast<operator_token*>(token_p)->operator_id()) {
					case ADD_OPERATOR: {
						if(infix_p_list.empty()) {
							delete token_p;
							token_p = nullptr;
							break;
						}
						if(token::check(infix_p_list.back(), OPEN_SYMBOL)) {
							delete token_p;
							token_p = nullptr;
							break;
						}
						if(token::check(infix_p_list.back(), ADD_OPERATOR)) {
							delete token_p;
							token_p = nullptr;
							break;
						}
						if(token::check(infix_p_list.back(), SUB_OPERATOR)) {
							delete token_p;
							token_p = nullptr;
							break;
						}
						if(token::check(infix_p_list.back(), NEG_OPERATOR)) {
							delete token_p;
							token_p = nullptr;
							break;
						}
						break;
					}
					case SUB_OPERATOR: {
						if(infix_p_list.empty()) {
							delete token_p;
							token_p = new operator_token(NEG_OPERATOR);
							break;
						}
						if(token::check(infix_p_list.back(), OPEN_SYMBOL)) {
							delete token_p;
							token_p = new operator_token(NEG_OPERATOR);
							break;
						}
						if(token::check(infix_p_list.back(), ADD_OPERATOR)) {
							delete infix_p_list.back();
							infix_p_list.back() = token_p;
							token_p = nullptr;
							break;
						}
						if(token::check(infix_p_list.back(), SUB_OPERATOR)) {
							delete infix_p_list.back();
							infix_p_list.back() =
								new operator_token(ADD_OPERATOR);
							delete token_p;
							token_p = nullptr;
							break;
						}
						if(token::check(infix_p_list.back(), NEG_OPERATOR)) {
							delete infix_p_list.back();
							infix_p_list.pop_back();
							delete token_p;
							token_p = nullptr;
							break;
						}
						break;
					}
					case MUL_OPERATOR: {
						if(infix_p_list.empty())
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"first token is mul operator"
							);
						if(token::check(infix_p_list.back(), OPERATOR_TOKEN))
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"consecutive operator tokens detected"
							);
						break;
					}
					case DIV_OPERATOR: {
						if(infix_p_list.empty())
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"first token is div operator"
							);
						if(token::check(infix_p_list.back(), OPERATOR_TOKEN))
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"consecutive operator tokens detected"
							);
						break;
					}
					case NEG_OPERATOR: {
						if(infix_p_list.empty())
							break;
						if(token::check(infix_p_list.back(), ADD_OPERATOR)) {
							delete token_p;
							token_p = nullptr;
							break;
						}
						if(token::check(infix_p_list.back(), SUB_OPERATOR)) {
							delete infix_p_list.back();
							infix_p_list.back() =
								new operator_token(ADD_OPERATOR);
							delete token_p;
							token_p = nullptr;
							break;
						}
						if(token::check(infix_p_list.back(), NEG_OPERATOR)) {
							delete infix_p_list.back();
							infix_p_list.pop_back();
							delete token_p;
							token_p = nullptr;
							break;
						}
						break;
					}
					case POW_OPERATOR: {
						if(infix_p_list.empty())
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"first token is pow operator"
							);
						if(token::check(infix_p_list.back(), OPERATOR_TOKEN))
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"consecutive operator tokens detected"
							);
						break;
					}
				}
				break;
			}
			case SYMBOL_TOKEN: {
				switch(dynamic_cast<symbol_token*>(token_p)->symbol_id()) {
					case OPEN_SYMBOL: {
						nesting_count++;
						if(infix_p_list.empty())
							break;
						if(token::check(infix_p_list.back(), CLOSE_SYMBOL))
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"open symbol after close symbol detected"
							);
						if(token::check(infix_p_list.back(), NUMERIC_TOKEN))
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"open symbol after numeric token detected"
							);
						if(token::check(infix_p_list.back(), VARIABLE_TOKEN))
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"open symbol after variable token detected"
							);
						break;
					}
					case CLOSE_SYMBOL: {
						nesting_count--;
						if(infix_p_list.empty())
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"first token is close symbol"
							);
						if(token::check(infix_p_list.back(), OPEN_SYMBOL))
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"close symbol after open symbol detected"
							);
						if(nesting_count < 0)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
								"parenthesis mismatch"
							);
						break;
					}
				}
				break;
			}
			case FUNCTION_TOKEN: {
				if(infix_p_list.empty())
					break;
				if(token::check(infix_p_list.back(), NUMERIC_TOKEN))
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"function operator after numeric token detected"
					);
				if(token::check(infix_p_list.back(), VARIABLE_TOKEN))
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"function operator after variable token detected"
					);
				if(token::check(infix_p_list.back(), CLOSE_SYMBOL))
					throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+
						"function operator after close symbol detected"
					);
				break;
			}
		}
		if(token_p != nullptr)
			infix_p_list.push_back(token_p);
	}
	if(infix_p_list.empty())
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ "empty token list");
	else if(nesting_count != 0)
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ "parenthesis mismatch");
	else if(token::check(infix_p_list.back(), OPERATOR_TOKEN))
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ "last token is an operator");
	else if(token::check(infix_p_list.back(), FUNCTION_TOKEN))
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ "last token is a function");
	return infix_to_postfix(infix_p_list);
}

std::list<mathexpr::token*> mathexpr::infix_to_postfix(std::list<token*>& infix_p_list) const {
	std::list<mathexpr::token*> postfix_p_list;
	std::stack<token*> token_p_stack;
	int nesting_count = 0;
	while(!infix_p_list.empty()) {
		token* token_p = infix_p_list.front();
		infix_p_list.pop_front();
		switch(token_p->token_id()) {
			case NUMERIC_TOKEN: {
				postfix_p_list.push_back(token_p);
				break;
			}
			case VARIABLE_TOKEN: {
				postfix_p_list.push_back(token_p);
				break;
			}
			case OPERATOR_TOKEN: {
				while(!token_p_stack.empty()) {
					if(token::check(token_p_stack.top(), OPEN_SYMBOL))
						break;
					const auto p1 = token_p->precedence();
					const auto p2 = token_p_stack.top()->precedence();
					if((p1 < 0) && ((std::fabs(p1) > std::fabs(p2))))
						break;
					if((p1 > 0) && ((std::fabs(p1) >= std::fabs(p2))))
						break;
					postfix_p_list.push_back(token_p_stack.top());
					token_p_stack.pop();
				}
				token_p_stack.push(token_p);
				break;
			}
			case SYMBOL_TOKEN: {
				switch(dynamic_cast<symbol_token*>(token_p)->symbol_id()) {
					case OPEN_SYMBOL: {
						token_p_stack.push(token_p);
						nesting_count++;
						break;
					}
					case CLOSE_SYMBOL: {
						nesting_count--;
						if(nesting_count < 0)
							throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ "parenthesis mismatch");
						while(!token::check(token_p_stack.top(), OPEN_SYMBOL)) {
							postfix_p_list.push_back(token_p_stack.top());
							token_p_stack.pop();
						}
						delete token_p_stack.top();
						token_p_stack.pop();
						delete token_p;
						break;
					}
				}
				break;
			}
			case FUNCTION_TOKEN: {
				token_p_stack.push(token_p);
				break;
			}
		}
	}
	while(!token_p_stack.empty()) {
		postfix_p_list.push_back(token_p_stack.top());
		token_p_stack.pop();
	}
	if(nesting_count != 0)
		throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+ "parenthesis mismatch");
	return postfix_p_list;
}

void mathexpr::delete_tokens() {
	for(auto* const token_p : _postfix_p_list)
		if(token_p != nullptr)
			delete token_p;
	_postfix_p_list.clear();
}

}