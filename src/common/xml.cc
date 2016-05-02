/**
 * @file src/common/xml.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <common/xml.hh>
#include <sstream>
#include <stack>
#include <cctype>

namespace xml {

element::element(std::istream& is) {
    const int _eof = std::char_traits<char>::eof();
    auto _peek_valid_char = [&]()->bool {
        const int c = is.peek();
        if(std::isalnum(c) || (c == '@') || (c == '-') || (c == '_'))
            return true;
        return false;
    };
    int row = 1;
    auto _pop_char = [&]()->int {
        if(is.peek() == '\n')
            row++;
        return is.get();
    };
    auto _pop_space = [&]() {
        while(std::isspace(is.peek()))
            _pop_char();
    };
    auto _throw_unexpected_eof = [&]() {
        throw std::runtime_error("unexpected end of file");
    };
    auto _throw_unexpected_char = [&](int c) {
        if(c == std::char_traits<char>::eof())
            _throw_unexpected_eof();
        std::ostringstream oss;
        oss << "unexpected character '" << static_cast<char>(c) << "'";
        oss << " in line " << row;
        throw std::runtime_error(oss.str());
    };
    auto _throw_unexpected_tag = [&](const std::string& tag) {
        std::ostringstream oss;
        oss << "unexpected tag '" << tag << "'";
        oss << " in line " << row;
        throw std::runtime_error(oss.str());
    };
    std::stack<element*> element_p_stack;
    element_p_stack.push(this);
    while(is.peek() != _eof) {
        while((is.peek() != _eof) && (is.peek() != '<'))
            _pop_char();
        _pop_char();
        _pop_space();
        if(is.peek() == '!') {
            /* comment */
            for(int i = 0; i < 2; i ++) {
                _pop_char();
                _pop_space();
                if(is.peek() != '-')
                    _throw_unexpected_char(is.peek());
            }
            for(int i = 0; i < 3; i++) {
                _pop_char();
                _pop_space();
                if(is.peek() == _eof)
                    _throw_unexpected_eof();
                else if(i < 2) {
                    if(is.peek() != '-')
                        i = -1;
                } else if(is.peek() != '>')
                    i = -1;
            }
        } else if(is.peek() == '/') {
            /* close tag */
            _pop_char();
            _pop_space();
            std::string name_str;
            while((is.peek() != _eof) && _peek_valid_char())
                name_str += _pop_char();
            if(name_str.empty())
                _throw_unexpected_char(is.peek());
            _pop_space();
            if(is.peek() != '>')
                _throw_unexpected_char(is.peek());
            _pop_char();
            if(element_p_stack.size() == 1)
                _throw_unexpected_tag(name_str);
            if(name_str != element_p_stack.top()->_name_str)
                _throw_unexpected_tag(name_str);
            element_p_stack.pop();
        } else if(is.peek() != _eof) {
            /* open or empty tag */
            std::string name_str;
            while((is.peek() != _eof) && _peek_valid_char())
                name_str += _pop_char();
            if(name_str.empty())
                _throw_unexpected_char(is.peek());
            _pop_space();
            element* child_p = new element;
            child_p->_name_str = name_str;
            child_p->_attr_map["xml::line"] = std::to_string(row);
            element_p_stack.top()->_child_p_vec.push_back(child_p);
            while((is.peek() != _eof) && (is.peek() != '/') && (is.peek() != '>')) {
                std::string attr_key;
                while((is.peek() != _eof) && _peek_valid_char())
                    attr_key += _pop_char();
                if(attr_key.empty())
                    _throw_unexpected_char(is.peek());
                _pop_space();
                if(is.peek() != '=')
                    _throw_unexpected_char(is.peek());
                _pop_char();
                _pop_space();
                if(is.peek() != '"')
                    _throw_unexpected_char(is.peek());
                _pop_char();
                std::string attr_value;
                while((is.peek() != _eof) && (is.peek() != '"'))
                    attr_value += _pop_char();
                if(is.peek() != '"')
                    _throw_unexpected_char(is.peek());
                _pop_char();
                _pop_space();
                child_p->_attr_map[attr_key] = attr_value;
            }
            if(is.peek() == '/') {
                _pop_char();
                _pop_space();
            } else
                element_p_stack.push(child_p);
            if(is.peek() != '>')
                _throw_unexpected_char(is.peek());
            _pop_char();
        }
    }
    if(element_p_stack.size() > 1)
        _throw_unexpected_eof();
}

element::element(const element& obj) {
    *this = obj;
}

element::~element() {
    for(const element* child_p : _child_p_vec)
        delete child_p;
    _child_p_vec.clear();
}

element& element::operator=(const element& obj) {
    if(this == &obj)
        return *this;
    for(const element* child_p : _child_p_vec)
        delete child_p;
    _child_p_vec.clear();
    _name_str = obj._name_str;
    _attr_map = obj._attr_map;
    for(const element* child_p : obj._child_p_vec)
        _child_p_vec.push_back(new element(*child_p));
    return *this;
}

int element::size() const {
    return _child_p_vec.size();
}

std::string element::name() const {
    return _name_str;
}

std::string element::attr(const std::string& key) const {
    auto cit = _attr_map.find(key);
    if(cit == _attr_map.cend())
        return std::string();
    return cit->second;
}

std::vector<const element*> element::children(const std::string& name_tag) const {
    std::vector<const element*> child_p_vec;
    for(const element* child_p : _child_p_vec)
        if(child_p->name() == name_tag)
            child_p_vec.push_back(child_p);
    return child_p_vec;
}


std::vector<const element*>::const_iterator element::cbegin() const {
    return _child_p_vec.cbegin();
}

std::vector<const element*>::const_iterator element::cend() const {
    return _child_p_vec.cend();
}

std::vector<const element*>::const_iterator element::begin() const {
    return _child_p_vec.cbegin();
}

std::vector<const element*>::const_iterator element::end() const {
    return _child_p_vec.cend();
}

};