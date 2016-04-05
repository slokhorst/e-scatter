/**
 * @file src/edge-detect/main.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <cfloat>
#include <iostream>
#include <stdexcept>
#include <string>
#include <boost/program_options.hpp>
#include <cpl/text.h>
#include "float_image.h"
#include "gauss2_model.h"
#include "nelder_mead.h"

namespace po = boost::program_options;

std::vector<double> fit_edges(const float_image& image, const std::function<double(const double&)>& model) {
	std::vector<double> edge_dx_vec(image.height());
	const nelder_mead::params<1> min_simplex = {-0.5};
	const nelder_mead::params<1> max_simplex = {+0.5};
	#pragma omp parallel for
	for(size_t k = 0; k < edge_dx_vec.size(); k++) {
		auto chi_2_fun = [&image, &model, &k](const nelder_mead::params<1>& params)->double {
			const double& dx = params[0];
			double chi_2 = 0;
			for(size_t i = 0; i < image.width(); i++)
				chi_2 += (image(i,k)-model(i-dx))*(image(i,k)-model(i-dx));
			return chi_2;
		};
		auto simplex = nelder_mead::randomize<1>(min_simplex, max_simplex);
		while(nelder_mead::iterate<1>(chi_2_fun, simplex) > 10*DBL_EPSILON);
		edge_dx_vec[k] = simplex[0][0];
	}
	return edge_dx_vec;
}

int main(int argc, char* argv[]) {
	srand48(time(nullptr));

	std::string input_image_file;
	bool transpose_image_flag = false;
	bool column_output_flag = false;

	std::string usage = "Usage: edge-detect [options] <input edge image>";
	po::options_description visible_opt("allowed options");
	visible_opt.add_options()
		("help,h", "produce help message")
		("transpose-image", "transpose the input image before processing")
		("column-output", "output the edge dispacements in a column instead of a row")
	;
	po::options_description hidden_opt("Hidden options");
	hidden_opt.add_options()
		("input-edge-image,i", po::value<std::string>(&input_image_file)->required(),
			"set input image file")
	;
	po::options_description options;
	options.add(visible_opt).add(hidden_opt);
	po::positional_options_description pos_opt;
	pos_opt.add("input-edge-image", 1);

	try {
		po::variables_map vm;
		po::store(po::command_line_parser(argc, argv).options(options).positional(pos_opt).run(), vm);
		if(vm.count("help")) {
			std::clog << usage << std::endl
			          << visible_opt << std::endl;
			return EXIT_SUCCESS;
		}
		po::notify(vm);

		if(vm.count("transpose-image"))
			transpose_image_flag = true;
		if(vm.count("column-output"))
			column_output_flag = true;
	} catch(std::exception& e) {
		std::clog << e.what() << std::endl;
		std::clog << usage << std::endl
		          << visible_opt << std::endl;
		return EXIT_FAILURE;
	}

	try {
		std::cout << cpl::text::cat("attempting to read image '", input_image_file, "'") << std::endl;
		float_image image(input_image_file);
		if(transpose_image_flag)
			image = image.transpose();
		std::cout << cpl::text::cat("image width = ", image.width(), " [px], height = ", image.height(), " [px]") << std::endl;
		std::cout << "starting profile fitting (Nelder-Mead simplex method)" << std::endl;
		const gauss2_model model(image);
		std::cout << cpl::text::cat("chi squared = ", cpl::text::float64(model.chi_2)) << std::endl;
		std::cout << cpl::text::cat(
			"profile center = ", cpl::text::float64(model.center), " [px]",
			", amplitude = ", cpl::text::float64(model.amplitude)
		) << std::endl;
		std::cout << cpl::text::cat(
			"profile first base = ", cpl::text::float64(model.base.first),
			", first std = ", cpl::text::float64(model.std.first), " [px]"
		) << std::endl;
		std::cout << cpl::text::cat(
			"profile second base = ", cpl::text::float64(model.base.second),
			", second std = ", cpl::text::float64(model.std.second), " [px]"
		) << std::endl;
		std::cout << "profile fitting completed" << std::endl;
		std::cout << "starting profile based edge detection (Nelder-Mead simplex method)" << std::endl;
		const auto edge_dx_vec = fit_edges(image, model);
		std::cout << "profile based edge detection completed" << std::endl;
		if(column_output_flag) {
			for(size_t i = 0; i < edge_dx_vec.size(); i++)
				std::cout << cpl::text::float64(edge_dx_vec[i]) << std::endl;
		} else {
			std::cout << cpl::text::float64(edge_dx_vec[0]);
			for(size_t i = 1; i < edge_dx_vec.size(); i++)
				std::cout << ", " << cpl::text::float64(edge_dx_vec[i]);
				std::cout << std::endl;
		}
		std::cout << std::flush;
		std::cout << "program exit" << std::endl;
		return EXIT_SUCCESS;
	} catch(const std::exception& e) {
		std::cout << "error:" << e.what() << std::endl;
		return EXIT_FAILURE;
	}
}