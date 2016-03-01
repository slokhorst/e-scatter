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
 * @file src/interpolate.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "interpolate.h"

double interpolate(const std::map<double,double>& xy_map, double x) {
	if(xy_map.empty())
		return 0;
	if(xy_map.size() == 1)
		return xy_map.cbegin()->second;
	std::map<double,double>::const_iterator cit;
	if(x <= xy_map.cbegin()->first)
		cit = std::next(xy_map.cbegin());
	else
		cit = xy_map.lower_bound(x);
	const double x2 = cit->first;
	const double y2 = cit->second;
	cit--;
	const double x1 = cit->first;
	const double y1 = cit->second;
	return y1+(y2-y1)*(x-x1)/(x2-x1);
}

double interpolate(const std::map<double,std::map<double,double>>& xyz_map, double x, double y) {
	if(xyz_map.empty())
		return 0;
	if(xyz_map.size() == 1)
		return interpolate(xyz_map.cbegin()->second, y);
	std::map<double,std::map<double,double>>::const_iterator cit;
	if(x <= xyz_map.cbegin()->first)
		cit = std::next(xyz_map.cbegin());
	else
		cit = xyz_map.lower_bound(x);
	const double x2 = cit->first;
	const double z2 = interpolate(cit->second, y);
	cit--;
	const double x1 = cit->first;
	const double z1 = interpolate(cit->second, y);
	return z1+(z2-z1)*(x-x1)/(x2-x1);
}