/**
 * @file src/line-gen/main.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
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
#include <cdsem/triangle.hh>
#include <cdsem/point3.hh>
#include <line-gen/vertex_grid.h>

namespace po = boost::program_options;

int8_t SUBS_MAT = 0, LINE_MAT = 0, vacuum_mat = -123;
int8_t term_mat = -127, mirror_mat = -122;
int8_t detector_se_mat = -125, detector_bs_mat = -124;

double PSD_Palasantzas2D(double s, double xi, double H, double k) {
	double p0 = 4*M_PI*s*s*xi*xi*H;
	return p0/pow(1+k*xi*k*xi,H+1);
}

std::vector<triangle> rect_to_tris(std::array<point3,4> verts) {
	return {
		triangle(verts[0],verts[1],verts[2],0,0),
		triangle(verts[0],verts[2],verts[3],0,0)
	};
}

std::vector<triangle> grid_to_tris(vertex_grid vg) {
	std::vector<triangle> tri_vec;
	for(int i=0; i<vg.u-2; i+=2) {
		for(int j=0; j<vg.v-2; j+=2) {
			tri_vec.push_back(triangle(vg.points[i+1][j+1],vg.points[i+2][j],vg.points[i][j],0,0));
			tri_vec.push_back(triangle(vg.points[i+1][j+1],vg.points[i+2][j+2],vg.points[i+2][j],0,0));
			tri_vec.push_back(triangle(vg.points[i+1][j+1],vg.points[i][j+2],vg.points[i+2][j+2],0,0));
			tri_vec.push_back(triangle(vg.points[i+1][j+1],vg.points[i][j],vg.points[i][j+2],0,0));
		}
	}
	return tri_vec;
}

std::string to_txt(const triangle& tri, int8_t imat, int8_t omat) {
	std::stringstream ss;
	ss.precision(6);
	ss << std::scientific;
	ss << int(imat) << ' ' << int(omat);
	ss << ' ' << tri.A.x*1e9<<' '<<tri.A.y*1e9<<' '<<tri.A.z*1e9;
	ss << ' ' << tri.B.x*1e9<<' '<<tri.B.y*1e9<<' '<<tri.B.z*1e9;
	ss << ' ' << tri.C.x*1e9<<' '<<tri.C.y*1e9<<' '<<tri.C.z*1e9;
	return ss.str();
}

class boundary {
public:
	boundary(std::vector<triangle> tris, int8_t imat, int8_t omat)
	: tri_vec(tris), in_mat(imat), out_mat(omat) {};
	std::vector<triangle> tri_vec;
	int8_t in_mat, out_mat;
	std::string to_txt() const {
		std::stringstream ss;
		for(const triangle& tri : tri_vec)
			ss << ::to_txt(tri, in_mat, out_mat) << std::endl;
		return ss.str();
	}
};

std::vector<boundary> make_line(
	double dx,
	double sx, double sy, double sz,
	double w, double h,
	vertex_grid vg_l, vertex_grid vg_r,
	point3 offset
) {
	auto vg_l_minmax = vg_l.get_z_minmax();
	vg_l_minmax.first = -w/2-(vg_l_minmax.first-1e-9);
	vg_l_minmax.second = -w/2-(vg_l_minmax.second+1e-9);
	auto vg_r_minmax = vg_r.get_z_minmax();
	vg_r_minmax.first = +w/2+(vg_r_minmax.first-1e-9);
	vg_r_minmax.second = +w/2+(vg_r_minmax.second+1e-9);
	std::vector<boundary> boundary_vec;
	//left wall
	vg_l.transform([w,h](point3 r) {
		return rotate(r,point3(0,-0.5*constant::pi,0)) + point3(-w/2,0,+h/2);
	});
	boundary_vec.push_back(boundary(grid_to_tris(vg_l), vacuum_mat, LINE_MAT));
	//right wall
	vg_r.transform([w,h](point3 r) {
		return rotate(r,point3(0,+0.5*constant::pi,0)) + point3(+w/2,0,+h/2);
	});
	boundary_vec.push_back(boundary(grid_to_tris(vg_r), vacuum_mat, LINE_MAT));
	//top
	{
		std::vector<triangle> tris;
		int i_l = vg_l.u-2;
		int i_r = 0;
		for(int j=0; j<vg_l.v-2; j+=2) {
			for(const triangle& tri : rect_to_tris({
				point3(vg_l_minmax.first,       vg_r.points[i_l][j+2].y, h),
				point3(vg_l_minmax.first,       vg_r.points[i_l][j].y,   h),
				point3(vg_l.points[i_l][j].x,   vg_l.points[i_l][j].y,   h),
				point3(vg_l.points[i_l][j+2].x, vg_l.points[i_l][j+2].y, h)
			}))
				tris.push_back(tri);
			for(const triangle& tri : rect_to_tris({
				point3(vg_r.points[i_r][j+2].x, vg_r.points[i_r][j+2].y, h),
				point3(vg_r.points[i_r][j].x,   vg_r.points[i_r][j].y,   h),
				point3(vg_r_minmax.first,       vg_l.points[i_r][j].y,   h),
				point3(vg_r_minmax.first,       vg_l.points[i_r][j+2].y, h)
			}))
				tris.push_back(tri);
		}
		for(const triangle& tri : rect_to_tris({
			point3(vg_r_minmax.first, +sy/2, h),
			point3(vg_r_minmax.first, -sy/2, h),
			point3(vg_l_minmax.first, -sy/2, h),
			point3(vg_l_minmax.first, +sy/2, h)
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, vacuum_mat, LINE_MAT));
	}
	//left substrate
	{
		std::vector<triangle> tris;
		int i=0;
		for(int j=0; j<vg_l.v-2; j+=2) {
			for(const triangle& tri : rect_to_tris({
				point3(vg_l.points[i][j+2].x, vg_l.points[i][j+2].y, 0),
				point3(vg_l.points[i][j].x,   vg_l.points[i][j].y,   0),
				point3(vg_l_minmax.second,    vg_l.points[i][j].y,   0),
				point3(vg_l_minmax.second,    vg_l.points[i][j+2].y, 0)
			}))
				tris.push_back(tri);
		}
		for(const triangle& tri : rect_to_tris({
			point3(-sx/2,              -sy/2, 0),
			point3(-sx/2,              +sy/2, 0),
			point3(vg_l_minmax.second, +sy/2, 0),
			point3(vg_l_minmax.second, -sy/2, 0)
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, vacuum_mat, SUBS_MAT));
	}
	//mid substrate
	{
		std::vector<triangle> tris;
		int i_l = 0;
		int i_r = vg_r.u-2;
		for(int j=0; j<vg_l.v-2; j+=2) {
			for(const triangle& tri : rect_to_tris({
				point3(vg_l_minmax.first,       vg_r.points[i_l][j+2].y, 0),
				point3(vg_l_minmax.first,       vg_r.points[i_l][j].y,   0),
				point3(vg_l.points[i_l][j].x,   vg_l.points[i_l][j].y,   0),
				point3(vg_l.points[i_l][j+2].x, vg_l.points[i_l][j+2].y, 0)
			}))
				tris.push_back(tri);
			for(const triangle& tri : rect_to_tris({
				point3(vg_r.points[i_r][j+2].x, vg_r.points[i_r][j+2].y, 0),
				point3(vg_r.points[i_r][j].x,   vg_r.points[i_r][j].y,   0),
				point3(vg_r_minmax.first,       vg_l.points[i_r][j].y,   0),
				point3(vg_r_minmax.first,       vg_l.points[i_r][j+2].y, 0)
			}))
				tris.push_back(tri);
		}
		for(const triangle& tri : rect_to_tris({
			point3(vg_r_minmax.first, +sy/2, 0),
			point3(vg_r_minmax.first, -sy/2, 0),
			point3(vg_l_minmax.first, -sy/2, 0),
			point3(vg_l_minmax.first, +sy/2, 0)
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, LINE_MAT, SUBS_MAT));
	}
	//right substrate
	{
		std::vector<triangle> tris;
		int i=vg_r.u-2;
		for(int j=0; j<vg_r.v-2; j+=2) {
			for(const triangle& tri : rect_to_tris({
				point3(vg_r.points[i][j].x,   vg_r.points[i][j].y,   0),
				point3(vg_r.points[i][j+2].x, vg_r.points[i][j+2].y, 0),
				point3(vg_r_minmax.second,    vg_r.points[i][j+2].y, 0),
				point3(vg_r_minmax.second,    vg_r.points[i][j].y,   0)
			}))
				tris.push_back(tri);
		}
		for(const triangle& tri : rect_to_tris({
			point3(+sx/2,              +sy/2, 0),
			point3(+sx/2,              -sy/2, 0),
			point3(vg_r_minmax.second, -sy/2, 0),
			point3(vg_r_minmax.second, +sy/2, 0)
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, vacuum_mat, SUBS_MAT));
	}

	// detectors
	{
		std::vector<triangle> se_tris;
		//-x SE
		for(const triangle& tri : rect_to_tris({
			point3(vg_l_minmax.second, -sy/2, +h+1e-9),
			point3(vg_l_minmax.second, -sy/2, 0      ),
			point3(vg_l_minmax.second, +sy/2, 0      ),
			point3(vg_l_minmax.second, +sy/2, +h+1e-9),
		}))
			se_tris.push_back(tri);
		//+x SE
		for(const triangle& tri : rect_to_tris({
			point3(vg_r_minmax.second, +sy/2, +h+1e-9),
			point3(vg_r_minmax.second, +sy/2, 0      ),
			point3(vg_r_minmax.second, -sy/2, 0      ),
			point3(vg_r_minmax.second, -sy/2, +h+1e-9),
		}))
			se_tris.push_back(tri);
		//+z SE
		for(const triangle& tri : rect_to_tris({
			point3(+sx/2, +sy/2, h+1e-9),
			point3(+sx/2, -sy/2, h+1e-9),
			point3(-sx/2, -sy/2, h+1e-9),
			point3(-sx/2, +sy/2, h+1e-9),
		}))
			se_tris.push_back(tri);
		boundary_vec.push_back(boundary(se_tris, detector_se_mat, vacuum_mat));
		//+z BSE
		std::vector<triangle> bs_tris;
		for(const triangle& tri : rect_to_tris({
			point3(+sx/2, +sy/2, h+2e-9),
			point3(+sx/2, -sy/2, h+2e-9),
			point3(-sx/2, -sy/2, h+2e-9),
			point3(-sx/2, +sy/2, h+2e-9),
		}))
			bs_tris.push_back(tri);
		boundary_vec.push_back(boundary(bs_tris, detector_bs_mat, vacuum_mat));
	}

	// kill planes
	{
		std::vector<triangle> tris;
		//+z
		for(const triangle& tri : rect_to_tris({
			point3(+sx/2, +sy/2, h+3e-9),
			point3(+sx/2, -sy/2, h+3e-9),
			point3(-sx/2, -sy/2, h+3e-9),
			point3(-sx/2, +sy/2, h+3e-9),
		}))
			tris.push_back(tri);
		//-z
		for(const triangle& tri : rect_to_tris({
			point3(+sx/2, -sy/2, -sz),
			point3(+sx/2, +sy/2, -sz),
			point3(-sx/2, +sy/2, -sz),
			point3(-sx/2, -sy/2, -sz),
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, term_mat, vacuum_mat));
	}

	// mirrors
	{
		std::vector<triangle> tris;
		//-x
		for(const triangle& tri : rect_to_tris({
			point3(-sx/2, -sy/2, +h+3e-9),
			point3(-sx/2, -sy/2, -sz    ),
			point3(-sx/2, +sy/2, -sz    ),
			point3(-sx/2, +sy/2, +h+3e-9),
		}))
			tris.push_back(tri);
		//+x
		for(const triangle& tri : rect_to_tris({
			point3(+sx/2, +sy/2, +h+3e-9),
			point3(+sx/2, +sy/2, -sz    ),
			point3(+sx/2, -sy/2, -sz    ),
			point3(+sx/2, -sy/2, +h+3e-9),
		}))
			tris.push_back(tri);
		//-y
		for(const triangle& tri : rect_to_tris({
			point3(+sx/2, -sy/2, +h+3e-9),
			point3(+sx/2, -sy/2, -sz    ),
			point3(-sx/2, -sy/2, -sz    ),
			point3(-sx/2, -sy/2, +h+3e-9),
		}))
			tris.push_back(tri);
		//+y
		for(const triangle& tri : rect_to_tris({
			point3(-sx/2, +sy/2, +h+3e-9),
			point3(-sx/2, +sy/2, -sz    ),
			point3(+sx/2, +sy/2, -sz    ),
			point3(+sx/2, +sy/2, +h+3e-9),
		}))
			tris.push_back(tri);
		boundary_vec.push_back(boundary(tris, mirror_mat, vacuum_mat));
	}

	for(boundary& b : boundary_vec)
		for(triangle& t : b.tri_vec) {
			t.A += offset;
			t.B += offset;
			t.C += offset;
		}

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
			std::clog
				<< "Generating roughness from Palasantzas' model with:\n"
				<< "\tσ = " << vars["std-dev"].as<double>() << " [m]\n"
				<< "\tξ = " << vars["cor-length"].as<double>() << " [m]\n"
				<< "\tH = " << vars["hurst-exp"].as<double>()
				<< std::endl;
		} else if(vars.count("import")) {
			method = IMPORT;
			import_fns = vars["import"].as<std::vector<std::string>>();
			if(import_fns.size() != 2*n_lines)
				throw std::runtime_error(
					"number of import walls ("+std::to_string(import_fns.size())+") "+
					"does not match the number of requested walls (2*"+std::to_string(n_lines)+")"
				);
			std::clog
				<< "Importing sidewall profile from:";
			for(const std::string& fn : import_fns)
				std::clog << "\n\t" << fn;
			std::clog << std::endl;
		} else {
			method = NONE;
			std::clog << "Generating smooth line." << std::endl;
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
			point3 offset(offset_x,0,0);

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
					if(!rgfs_l.is_open()) throw std::runtime_error("failed to open file '"+import_fns[2*i+0]+"'");
					if(!rgfs_r.is_open()) throw std::runtime_error("failed to open file '"+import_fns[2*i+1]+"'");
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
				std::cout << s.to_txt() << std::endl;

		}
	} catch(const std::exception& e) {
		std::clog << "Error applying roughness: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
