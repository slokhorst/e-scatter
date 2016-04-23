/**
 * @file src/sample-viewer/main.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <iostream>
#include <sstream>
#include <vector>
#include <sample-viewer/sample_viewer.h>
#include <cdsem/point3.hh>
#include <cdsem/triangle.hh>

int main(int argc, char* argv[]) {
	std::vector<triangle> mi_vec;

	uint ln=0;
	std::string line;
	try {
		while (std::getline(std::cin, line)) {
			ln++;
			if(line[0] == '#' || line.size() == 0)
				continue;
			std::stringstream ss(line);

			int in, out;
			double Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz;

			ss >> in >> out >> Ax >> Ay >> Az >> Bx >> By >> Bz >> Cx >> Cy >> Cz;
			mi_vec.push_back(
				triangle(point3(Ax, Ay, Az), point3(Bx, By, Bz), point3(Cx, Cy, Cz), in, out)
			);
		}
	} catch(std::exception& e) {
		std::clog << "error parsing geometry on line " << ln << ": " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	sample_viewer sv(&mi_vec);

	return EXIT_SUCCESS;
}
