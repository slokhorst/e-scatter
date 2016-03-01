/*
 * Copyright 2014-2016 Thomas Verduin <T.Verduin@tudelft.nl>
 *
 * This program is part of the electron-matter interaction program (SCATTER).
 *
 * SCATTER is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

/**
 * @file src/xml.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SCATTER__XML__HEADER_INCLUDED
#define SCATTER__XML__HEADER_INCLUDED

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
