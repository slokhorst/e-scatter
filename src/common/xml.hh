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
#include <memory>

namespace xml {

    class element {
    public:
        using value_type = const element *;
        using const_iterator = std::vector<value_type>::const_iterator;

        element(std::istream&);
        element(const element&);
        ~element();

        element& operator=(const element&);
        int size() const;
        std::string name() const;
        std::string attr(const std::string&) const;

        std::vector<value_type> children(const std::string& name_tag) const;
        const_iterator cbegin() const;
        const_iterator cend() const;
        const_iterator begin() const;
        const_iterator end() const;

    private:
        element() = default;
        std::string _name_str;
        std::map<std::string,std::string> _attr_map;
        std::vector<value_type> _child_p_vec;
    };

}

#endif
