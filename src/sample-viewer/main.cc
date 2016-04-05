/**
 * @file src/sample-viewer/main.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <sample-viewer/sample_viewer.h>
#include <iostream>
#include <sstream>

int main(int argc, char* argv[]) {
	std::vector<material_interface> mi_vec;

	uint ln=0;
	std::string line;
	try {
		while (std::getline(std::cin, line)) {
			ln++;
			if(line[0] == '#' || line.size() == 0)
				continue;
			std::stringstream ss(line);

			int in_mat, out_mat;
			double x1, y1, z1, x2, y2, z2, x3, y3, z3;

			ss >> in_mat >> out_mat >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
			mi_vec.push_back(
				material_interface(
					cpl::triangle({
						cpl::vector3(x1, y1, z1),
						cpl::vector3(x2, y2, z2),
						cpl::vector3(x3, y3, z3)
					}), in_mat, out_mat
				)
			);
		}
	} catch(std::exception& e) {
		std::clog << "error parsing geometry on line " << ln << ": " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	sample_viewer sv(&mi_vec);

	return EXIT_SUCCESS;
}
