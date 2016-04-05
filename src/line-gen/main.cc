/*
 * src/line-gen/main.cc
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

#include <algorithm>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <utility>
#include <boost/program_options.hpp>
#include <common/constant.hh>
#include <cpl/text.h>
#include <cpl/triangle.h>
#include <cpl/vector3.h>

#include <line-gen/vertex_grid.h>

namespace po = boost::program_options;

int8_t SUBS_MAT = 0, LINE_MAT = 0, vacuum_mat = -123;
int8_t term_mat = -127, mirror_mat = -122;
int8_t detector_se_mat = -125, detector_bs_mat = -124;

double PSD_Palasantzas2D(double s, double xi, double H, double k) {
	double p0 = 4*M_PI*s*s*xi*xi*H;
	return p0/pow(1+k*xi*k*xi,H+1);
}

std::vector<cpl::triangle> rect_to_tris(std::array<cpl::vector3,4> verts) {
	return {
		cpl::triangle({verts[0],verts[1],verts[2]}),
		cpl::triangle({verts[0],verts[2],verts[3]})
	};
}

std::vector<cpl::triangle> grid_to_tris(vertex_grid vg) {
	std::vector<cpl::triangle> tri_vec;
	for(int i=0; i<vg.u-2; i+=2) {
		for(int j=0; j<vg.v-2; j+=2) {
			tri_vec.push_back(cpl::triangle({vg.points[i+1][j+1],vg.points[i+2][j],vg.points[i][j]}));
			tri_vec.push_back(cpl::triangle({vg.points[i+1][j+1],vg.points[i+2][j+2],vg.points[i+2][j]}));
			tri_vec.push_back(cpl::triangle({vg.points[i+1][j+1],vg.points[i][j+2],vg.points[i+2][j+2]}));
			tri_vec.push_back(cpl::triangle({vg.points[i+1][j+1],vg.points[i][j],vg.points[i][j+2]}));
		}
	}
	return tri_vec;
}

std::string to_txt(const cpl::triangle& tri, int8_t imat, int8_t omat, double dz) {
	std::stringstream ss;
	ss.precision(6);
	ss << std::scientific;
	ss << int(imat) << ' ' << int(omat);
	ss << ' ' << tri.vertices()[0].x*1e9<<' '<<tri.vertices()[0].y*1e9<<' '<<(tri.vertices()[0].z+dz)*1e9;
	ss << ' ' << tri.vertices()[1].x*1e9<<' '<<tri.vertices()[1].y*1e9<<' '<<(tri.vertices()[1].z+dz)*1e9;
	ss << ' ' << tri.vertices()[2].x*1e9<<' '<<tri.vertices()[2].y*1e9<<' '<<(tri.vertices()[2].z+dz)*1e9;
	return ss.str();
}

class boundary {
public:
	boundary(std::vector<cpl::triangle> tris, int8_t imat, int8_t omat)
	: tri_vec(tris), in_mat(imat), out_mat(omat) {};
	std::vector<cpl::triangle> tri_vec;
	int8_t in_mat, out_mat;
	std::string to_txt(double dz) const {
		std::stringstream ss;
		for(const cpl::triangle& tri : tri_vec)
			ss << ::to_txt(tri, in_mat, out_mat, dz) << std::endl;
		return ss.str();
	}
};

std::vector<boundary> make_line(
	double dx,
	double sx, double sy, double sz,
	double w, double h,
	vertex_grid vg_l, vertex_grid vg_r,
	cpl::vector3 offset
) {
	auto vg_l_minmax = vg_l.get_z_minmax();
	vg_l_minmax.first = -w/2-(vg_l_minmax.first-1e-9);
	vg_l_minmax.second = -w/2-(vg_l_minmax.second+1e-9);
	auto vg_r_minmax = vg_r.get_z_minmax();
	vg_r_minmax.first = +w/2+(vg_r_minmax.first-1e-9);
	vg_r_minmax.second = +w/2+(vg_r_minmax.second+1e-9);
	std::vector<boundary> boundary_vec;
	//left wall
	vg_l.transform([w,h](cpl::vector3 r) {
		return r.rotate(cpl::vector3(0,-0.5*constant::pi,0)) + cpl::vector3(-w/2,0,+h/2);
	});
	boundary_vec.push_back(boundary(grid_to_tris(vg_l), vacuum_mat, LINE_MAT));
	//right wall
	vg_r.transform([w,h](cpl::vector3 r) {
		return r.rotate(cpl::vector3(0,+0.5*constant::pi,0)) + cpl::vector3(+w/2,0,+h/2);
	});
	boundary_vec.push_back(boundary(grid_to_tris(vg_r), vacuum_mat, LINE_MAT));
	//top
	{
		std::vector<cpl::triangle> tris;
		int i_l = vg_l.u-2;
		int i_r = 0;
		for(int j=0; j<vg_l.v-2; j+=2) {
			for(const cpl::triangle& tri : rect_to_tris({
				cpl::vector3(vg_l_minmax.first,       vg_r.points[i_l][j+2].y, h),
				cpl::vector3(vg_l_minmax.first,       vg_r.points[i_l][j].y,   h),
				cpl::vector3(vg_l.points[i_l][j].x,   vg_l.points[i_l][j].y,   h),
				cpl::vector3(vg_l.points[i_l][j+2].x, vg_l.points[i_l][j+2].y, h)
			}))
				tris.push_back(tri);
			for(const cpl::triangle& tri : rect_to_tris({
				cpl::vector3(vg_r.points[i_r][j+2].x, vg_r.points[i_r][j+2].y, h),
				cpl::vector3(vg_r.points[i_r][j].x,   vg_r.points[i_r][j].y,   h),
				cpl::vector3(vg_r_minmax.first,       vg_l.points[i_r][j].y,   h),
				cpl::vector3(vg_r_minmax.first,       vg_l.points[i_r][j+2].y, h)
			}))
				tris.push_back(tri);
		}
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(vg_r_minmax.first, +sy/2, h),
			cpl::vector3(vg_r_minmax.first, -sy/2, h),
			cpl::vector3(vg_l_minmax.first, -sy/2, h),
			cpl::vector3(vg_l_minmax.first, +sy/2, h)
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, vacuum_mat, LINE_MAT));
	}
	//left substrate
	{
		std::vector<cpl::triangle> tris;
		int i=0;
		for(int j=0; j<vg_l.v-2; j+=2) {
			for(const cpl::triangle& tri : rect_to_tris({
				cpl::vector3(vg_l.points[i][j+2].x, vg_l.points[i][j+2].y, 0),
				cpl::vector3(vg_l.points[i][j].x,   vg_l.points[i][j].y,   0),
				cpl::vector3(vg_l_minmax.second,    vg_l.points[i][j].y,   0),
				cpl::vector3(vg_l_minmax.second,    vg_l.points[i][j+2].y, 0)
			}))
				tris.push_back(tri);
		}
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(-sx/2,              -sy/2, 0),
			cpl::vector3(-sx/2,              +sy/2, 0),
			cpl::vector3(vg_l_minmax.second, +sy/2, 0),
			cpl::vector3(vg_l_minmax.second, -sy/2, 0)
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, vacuum_mat, SUBS_MAT));
	}
	//mid substrate
	{
		std::vector<cpl::triangle> tris;
		int i_l = 0;
		int i_r = vg_r.u-2;
		for(int j=0; j<vg_l.v-2; j+=2) {
			for(const cpl::triangle& tri : rect_to_tris({
				cpl::vector3(vg_l_minmax.first,       vg_r.points[i_l][j+2].y, 0),
				cpl::vector3(vg_l_minmax.first,       vg_r.points[i_l][j].y,   0),
				cpl::vector3(vg_l.points[i_l][j].x,   vg_l.points[i_l][j].y,   0),
				cpl::vector3(vg_l.points[i_l][j+2].x, vg_l.points[i_l][j+2].y, 0)
			}))
				tris.push_back(tri);
			for(const cpl::triangle& tri : rect_to_tris({
				cpl::vector3(vg_r.points[i_r][j+2].x, vg_r.points[i_r][j+2].y, 0),
				cpl::vector3(vg_r.points[i_r][j].x,   vg_r.points[i_r][j].y,   0),
				cpl::vector3(vg_r_minmax.first,       vg_l.points[i_r][j].y,   0),
				cpl::vector3(vg_r_minmax.first,       vg_l.points[i_r][j+2].y, 0)
			}))
				tris.push_back(tri);
		}
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(vg_r_minmax.first, +sy/2, 0),
			cpl::vector3(vg_r_minmax.first, -sy/2, 0),
			cpl::vector3(vg_l_minmax.first, -sy/2, 0),
			cpl::vector3(vg_l_minmax.first, +sy/2, 0)
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, LINE_MAT, SUBS_MAT));
	}
	//right substrate
	{
		std::vector<cpl::triangle> tris;
		int i=vg_r.u-2;
		for(int j=0; j<vg_r.v-2; j+=2) {
			for(const cpl::triangle& tri : rect_to_tris({
				cpl::vector3(vg_r.points[i][j].x,   vg_r.points[i][j].y,   0),
				cpl::vector3(vg_r.points[i][j+2].x, vg_r.points[i][j+2].y, 0),
				cpl::vector3(vg_r_minmax.second,    vg_r.points[i][j+2].y, 0),
				cpl::vector3(vg_r_minmax.second,    vg_r.points[i][j].y,   0)
			}))
				tris.push_back(tri);
		}
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(+sx/2,              +sy/2, 0),
			cpl::vector3(+sx/2,              -sy/2, 0),
			cpl::vector3(vg_r_minmax.second, -sy/2, 0),
			cpl::vector3(vg_r_minmax.second, +sy/2, 0)
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, vacuum_mat, SUBS_MAT));
	}

	// detectors
	{
		std::vector<cpl::triangle> se_tris;
		//-x SE
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(vg_l_minmax.second, -sy/2, +h+1e-9),
			cpl::vector3(vg_l_minmax.second, -sy/2, 0      ),
			cpl::vector3(vg_l_minmax.second, +sy/2, 0      ),
			cpl::vector3(vg_l_minmax.second, +sy/2, +h+1e-9),
		}))
			se_tris.push_back(tri);
		//+x SE
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(vg_r_minmax.second, +sy/2, +h+1e-9),
			cpl::vector3(vg_r_minmax.second, +sy/2, 0      ),
			cpl::vector3(vg_r_minmax.second, -sy/2, 0      ),
			cpl::vector3(vg_r_minmax.second, -sy/2, +h+1e-9),
		}))
			se_tris.push_back(tri);
		//+z SE
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(+sx/2, +sy/2, h+1e-9),
			cpl::vector3(+sx/2, -sy/2, h+1e-9),
			cpl::vector3(-sx/2, -sy/2, h+1e-9),
			cpl::vector3(-sx/2, +sy/2, h+1e-9),
		}))
			se_tris.push_back(tri);
		boundary_vec.push_back(boundary(se_tris, detector_se_mat, vacuum_mat));
		//+z BSE
		std::vector<cpl::triangle> bs_tris;
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(+sx/2, +sy/2, h+2e-9),
			cpl::vector3(+sx/2, -sy/2, h+2e-9),
			cpl::vector3(-sx/2, -sy/2, h+2e-9),
			cpl::vector3(-sx/2, +sy/2, h+2e-9),
		}))
			bs_tris.push_back(tri);
		boundary_vec.push_back(boundary(bs_tris, detector_bs_mat, vacuum_mat));
	}

	// kill planes
	{
		std::vector<cpl::triangle> tris;
		//+z
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(+sx/2, +sy/2, h+3e-9),
			cpl::vector3(+sx/2, -sy/2, h+3e-9),
			cpl::vector3(-sx/2, -sy/2, h+3e-9),
			cpl::vector3(-sx/2, +sy/2, h+3e-9),
		}))
			tris.push_back(tri);
		//-z
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(+sx/2, -sy/2, -sz),
			cpl::vector3(+sx/2, +sy/2, -sz),
			cpl::vector3(-sx/2, +sy/2, -sz),
			cpl::vector3(-sx/2, -sy/2, -sz),
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, term_mat, vacuum_mat));
	}

	// mirrors
	{
		std::vector<cpl::triangle> tris;
		//-x
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(-sx/2, -sy/2, +h+3e-9),
			cpl::vector3(-sx/2, -sy/2, -sz    ),
			cpl::vector3(-sx/2, +sy/2, -sz    ),
			cpl::vector3(-sx/2, +sy/2, +h+3e-9),
		}))
			tris.push_back(tri);
		//+x
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(+sx/2, +sy/2, +h+3e-9),
			cpl::vector3(+sx/2, +sy/2, -sz    ),
			cpl::vector3(+sx/2, -sy/2, -sz    ),
			cpl::vector3(+sx/2, -sy/2, +h+3e-9),
		}))
			tris.push_back(tri);
		//-y
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(+sx/2, -sy/2, +h+3e-9),
			cpl::vector3(+sx/2, -sy/2, -sz    ),
			cpl::vector3(-sx/2, -sy/2, -sz    ),
			cpl::vector3(-sx/2, -sy/2, +h+3e-9),
		}))
			tris.push_back(tri);
		//+y
		for(const cpl::triangle& tri : rect_to_tris({
			cpl::vector3(-sx/2, +sy/2, +h+3e-9),
			cpl::vector3(-sx/2, +sy/2, -sz    ),
			cpl::vector3(+sx/2, +sy/2, -sz    ),
			cpl::vector3(+sx/2, +sy/2, +h+3e-9),
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, mirror_mat, vacuum_mat));
	}

	for(boundary& b : boundary_vec)
		for(cpl::triangle& t : b.tri_vec)
			for(cpl::vector3& v : t.vertices())
				v += offset;

	return boundary_vec;
}

int main(int argc, char* argv[]) {
	double dx;
	double substrate_sx;
	double substrate_sy;
	double substrate_sz;
	double line_width, line_height;
	uint n_lines;

	enum roughness_method { NONE, IMPORT, THORSOS };
	enum roughness_method method = NONE;

	std::vector<std::string> import_fns;
	double s, xi, H;

	std::string usage = "Usage: line-gen [options]";
	po::options_description options("allowed options");
	options.add_options()
		("help,h", "produce help message")
		("dx", po::value<double>(&dx)->default_value(1e-9,"1e-9"),
			"set (major) grid spacing")
		("line-width", po::value<double>(&line_width)->default_value(32e-9,"32e-9"),
			"set line width")
		("line-height", po::value<double>(&line_height)->default_value(32e-9,"32e-9"),
			"set line height")
		("substrate-length", po::value<double>(&substrate_sy)->default_value(1000e-9,"1000e-9"),
			"set line length")
		("n-lines", po::value<uint>(&n_lines)->default_value(1),
			"set number of adjacent lines")

		("std-dev,s", po::value<double>(&s)->default_value(1e-9,"1e-9"),
			"set standard deviation σ for Palasantzas")
		("cor-length,x", po::value<double>(&xi)->default_value(20e-9,"20e-9"),
			"set correlation length ξ for Palasantzas")
		("hurst-exp,H", po::value<double>(&H)->default_value(0.75),
			"set Hurst exponent H for Palasantzas")

		("import", po::value<std::vector<std::string>>(),
			"set input roughness (use this argument once for every wall)")
	;

	try {
		po::variables_map vars;
		po::store(po::command_line_parser(argc, argv).options(options).run(), vars);
		if(vars.count("help")) {
			std::clog << usage << std::endl
			          << options << std::endl;
			return EXIT_SUCCESS;
		}
		po::notify(vars);

		if(!vars["std-dev"].defaulted() || !vars["cor-length"].defaulted() || !vars["hurst-exp"].defaulted()) {
			method = THORSOS;
			std::clog << cpl::text::cat(
				"Generating roughness from Palasantzas' model with:\n",
				"\t","σ = ",vars["std-dev"].as<double>()," [m]\n",
				"\t","ξ = ",vars["cor-length"].as<double>()," [m]\n",
				"\t","H = ",vars["hurst-exp"].as<double>()
			) << std::endl;
		} else if(vars.count("import")) {
			method = IMPORT;
			import_fns = vars["import"].as<std::vector<std::string>>();
			if(import_fns.size() != 2*n_lines)
				throw std::runtime_error(cpl::text::cat(
					"number of import walls (",import_fns.size(),") ",
					"does not match the number of requested walls (2*",n_lines,")"
				));
			std::clog << cpl::text::cat(
				"Importing sidewall profile from:");
			for(const std::string& fn : import_fns)
				std::clog << "\n\t" << fn;
			std::clog << std::endl;
		} else {
			method = NONE;
			std::clog << cpl::text::cat("Generating smooth line.") << std::endl;
		}
	} catch(const std::exception& e) {
		std::clog << e.what() << std::endl;
		std::clog << usage << std::endl
		          << options << std::endl;
		return EXIT_FAILURE;
	}

	substrate_sx = 2*line_width;
	substrate_sz = 1e-6;

	try {
		for (uint i=0; i<n_lines; i++) {
			double offset_x = (2*line_width) * (i - ((double)n_lines-1)/2.0);
			cpl::vector3 offset(offset_x,0,0);

			vertex_grid roughness_l(std::round(line_height/dx+1)*2, std::round(substrate_sy/dx+1)*2, dx/2, dx/2);
			vertex_grid roughness_r(std::round(line_height/dx+1)*2, std::round(substrate_sy/dx+1)*2, dx/2, dx/2);

			switch (method) {
				case THORSOS:
				{
					std::function<double(double)> PSD = [s,xi,H](double k) { return PSD_Palasantzas2D(s,xi,H,k); };
					roughness_l.set_z_thorsos(PSD, s);
					roughness_r.set_z_thorsos(PSD, s);
				}
				break;
				case IMPORT:
				{
					std::ifstream rgfs_l(import_fns[2*i+0]);
					std::ifstream rgfs_r(import_fns[2*i+1]);
					if(!rgfs_l.is_open()) throw std::runtime_error(cpl::text::cat("failed to open file '",import_fns[2*i+0],"'"));
					if(!rgfs_r.is_open()) throw std::runtime_error(cpl::text::cat("failed to open file '",import_fns[2*i+1],"'"));
					roughness_l.set_z_csv(rgfs_l);
					roughness_r.set_z_csv(rgfs_r);
				}
				break;
				case NONE:
				break;
			}

			std::vector<boundary> boundary_vec = make_line(
				dx,
				substrate_sx, substrate_sy, substrate_sz,
				line_width, line_height,
				roughness_l, roughness_r,
				offset
			);
			for(const boundary& s : boundary_vec)
				std::cout << s.to_txt(line_height-32e-9) << std::endl;

		}
	} catch(const std::exception& e) {
		std::clog << cpl::text::cat("Error applying roughness: ",e.what()) << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
