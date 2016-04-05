/*
 * src/line-gen/vertex_grid.h
 *
 * Copyright 2015 Thomas Verduin <T.Verduin@tudelft.nl>
 *                Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 *
 * This program is free software; you can redistribute it and/or modify
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
 *
 *
 */

#ifndef LINE_GEN__GRID__HEADER_INCLUDED
#define LINE_GEN__GRID__HEADER_INCLUDED

#include <iostream>
#include <functional>
#include <vector>

#include <cpl/vector3.h>

class vertex_grid
{
public:
	std::vector<std::vector<cpl::vector3>> points;
	int u, v;
	double dx, dy;

	vertex_grid() = default;
	vertex_grid(int _u, int _v, double u_spacing, double v_spacing);

	void transform(std::function<cpl::vector3(cpl::vector3)> transform);

	void set_z_csv(std::istream& input);
	void set_z_thorsos(std::function<double(double)> PSD, double sigma=-1);

	std::pair<double,double> get_z_minmax() const;

	void save_gnusurf(std::ostream& output) const;
	void save_matlab(std::ostream& output) const;
};

#endif