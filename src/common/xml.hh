/**
 * @file src/common/xml.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__XML__HEADER_INCLUDED
#define eSCATTER__COMMON__XML__HEADER_INCLUDED

#include <map>
#include <string>
#include <vector>

namespace xml {

class element {
public:
    element(std::istream&);
    element(const element&);
    ~element();
    element& operator=(const element&);
    int size() const;
    std::string name() const;
    std::string attr(const std::string&) const;
    std::vector<const element*>::const_iterator cbegin() const;
    std::vector<const element*>::const_iterator cend() const;
    std::vector<const element*>::const_iterator begin() const;
    std::vector<const element*>::const_iterator end() const;
private:
    element() = default;
    std::string _name_str;
    std::map<std::string,std::string> _attr_map;
    std::vector<const element*> _child_p_vec;
};

};

#endif