/*
 * Copyright 2014-20156 Thomas Verduin <T.Verduin@tudelft.nl>
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
 * @file src/interpolate.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef SCATTER__INTERPOLATE__HEADER_INCLUDED
#define SCATTER__INTERPOLATE__HEADER_INCLUDED

#include <map>

/*!
 * Linear interpolation in one dimension.
 * Computational complexity is O(log N).
 */
double interpolate(const std::map<double,double>& xy_map, double x);
/*!
 * Linear interpolation in two dimensions.
 * Computational complexity is O(log N).
 */
double interpolate(const std::map<double,std::map<double,double>>& xyz_map, double x, double y);

#endif