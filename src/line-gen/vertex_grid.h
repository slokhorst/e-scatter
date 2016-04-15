/**
 * @file src/line-gen/vertex-grid.h
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
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