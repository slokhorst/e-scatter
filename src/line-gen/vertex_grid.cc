/**
 * @file src/line-gen/vertex-grid.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <line-gen/vertex_grid.h>
#include <complex>
#include <random>
#include <cinttypes>
#include <fftw3.h>

std::vector<std::string> split_string(const std::string& str, char c) {
	std::vector<std::string> str_vec;
	std::string sub_str;
	for(auto cit = str.cbegin(); cit < str.cend(); cit++)
		if(*cit != c)
			sub_str.push_back(*cit);
		else {
			str_vec.push_back(sub_str);
			sub_str.clear();
		}
	if(!sub_str.empty())
		str_vec.push_back(sub_str);
	return str_vec;
}

vertex_grid::vertex_grid(int _u, int _v, double u_spacing, double v_spacing) {
	u = _u;
	v = _v;
	dx = u_spacing;
	dy = v_spacing;

	const double x0 = -(u-1)/2 * dx;
	const double y0 = -(v-1)/2 * dy;
	for (int i=0; i<u; i++) {
		points.push_back(std::vector<point3>(v));
		for (int j=0; j<v; j++) {
			points[i][j] = point3(x0 + i*dx, y0 + j*dy, 0);
		}
	}
}

void vertex_grid::transform(std::function<point3(point3)> transfunc) {
	for(int i=0; i<u; i++)
		for(int j=0; j<v; j++)
			points[i][j] = transfunc(points[i][j]);
}

void vertex_grid::set_z_csv(std::istream& input) {
	int i = 0;
	while(!input.eof())
	{
		std::string str;
		std::getline(input, str);
		//str = cpl::text::strip_string(str);
		if (str.empty())
			continue;

		const auto substr_vec = split_string(str, ',');
		if (substr_vec.size() == 0)
			continue;

		if (i >= u)
			throw std::runtime_error("cannot load input with u>"+std::to_string(i)+" into grid with u="+std::to_string(u));
		if (substr_vec.size() > (uint)v)
			throw std::runtime_error("cannot load input with v="+std::to_string(substr_vec.size())+" into grid with v="+std::to_string(v));

		for (uint j=0; j<substr_vec.size(); j++) {
			//TODO: fix
			//points[i][j].z = substr_vec[j];
		}
		i++;
	}
}

void vertex_grid::set_z_thorsos(std::function<double(double)> PSD, double sigma) {
	if(u%2!=0 || v%2!=0)
		throw std::runtime_error("u and v must be even for Thorsos");

	fftw_complex* data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * u * v);
	fftw_plan p = fftw_plan_dft_2d(u, v, data, data, FFTW_BACKWARD, FFTW_ESTIMATE);

	std::random_device rnd_dev;
	std::mt19937_64 rng(rnd_dev());
	std::normal_distribution<double> gaussian(0, 1);

	double Lx = dx*(u-1);
	double Ly = dy*(v-1);

	std::complex<double>** F = new std::complex<double>*[u];
	for(int i=0; i<u; ++i)
		F[i] = new std::complex<double>[v];

	for(int i=(-u/2); i<=(u/2-1); i++) {
	for(int j=(-v/2); j<=(v/2-1); j++) {
		double kx = (2*M_PI)/Lx * i;
		double ky = (2*M_PI)/Ly * j;
		std::complex<double> Y = std::complex<double>(gaussian(rng),gaussian(rng))/sqrt(2);
		if((i==0 && j==0) || i==(-u/2) || j==(-v/2))
			Y = std::complex<double>(gaussian(rng), 0);
		F[u/2+i][v/2+j] = Y*sqrt(Lx*Ly*PSD(sqrt(kx*kx+ky*ky)));
	}
	}

	for(int i=0; i<=(u/2-1); i++) {
	for(int j=0; j<=(v/2-1); j++) {
		F[u/2+i][v/2+j] = std::conj(F[u/2-i][v/2-j]);
		F[u/2-i][v/2+j] = std::conj(F[u/2+i][v/2-j]);
	}
	}
	for(int i=0; i<=(u/2-1); i++)
		F[u/2+i][0] = std::conj(F[u/2-i][0]);
	for(int j=0; j<=(v/2-1); j++)
		F[0][v/2+j] = std::conj(F[0][v/2-j]);


	for(int i=0; i<u; i++) {
	for(int j=0; j<v; j++) {
		// flip the 1st & 3rd, and the 2nd & 4th quadrant,
		// so the 0-frequency at (i,j)=(u/2,v/2) is at (is,js)=(0,0)
		int is = (i+u/2)%u;
		int js = (j+v/2)%v;

		data[v*is+js][0] = std::real(F[i][j]);
		data[v*is+js][1] = std::imag(F[i][j]);
	}
	}

	for(int i=0; i<u; ++i)
		delete[] F[i];
	delete[] F;

	fftw_execute(p);

	double sum = 0;
	for(int j=0; j<v; j++) {
	for(int i=0; i<u; i++) {
		points[i][j].z = data[v*i+j][0]/(Lx*Ly);
		sum += points[i][j].z;
	}
	}
	double mean = sum/(u*v);
	double sum_var = 0;
	for(int j=0; j<v; j++) {
	for(int i=0; i<u; i++) {
		points[i][j].z -= mean;
		sum_var += points[i][j].z*points[i][j].z;
	}
	}
	double var = sum_var/(u*v);
	double sigma_thorsos = sqrt(var);

	fftw_destroy_plan(p);
	fftw_free(data);

	if (sigma > 0) {
		// scale z
		// note: mean(z) == 0
		//       std(z) == sigma
		//       std(std(z)) == 0
		//   but std( xi(z)) != 0
		//   and std( rx(z)) != 0
		for(int i=0; i<u; i++) {
		for(int j=0; j<v; j++) {
			points[i][j].z *= sigma/sigma_thorsos;
		}
		}
	}
}

void vertex_grid::save_gnusurf(std::ostream& output) const {
	for (int i=0; i<u; i++) {
	for (int j=0; j<v; j++) {
		output
			<< (points[i][j].x) << ' '
			<< (points[i][j].y) << ' '
			<< (points[i][j].z) << ' '
			<< std::endl;
	}
		output << std::endl;
	}
}

std::pair<double,double> vertex_grid::get_z_minmax() const {
	std::pair<double,double> minmax(0,0);
	for (int i=0; i<u; i++) {
	for (int j=0; j<v; j++) {
		minmax.first = std::min(minmax.first,points[i][j].z);
		minmax.second = std::max(minmax.second,points[i][j].z);
	}
	}
	return minmax;
}


void vertex_grid::save_matlab(std::ostream& output) const {
	for (int i=0; i<u; i++) {
	for (int j=0; j<v; j++) {
		output << (points[i][j].z) << ' ';
	}
		output << std::endl;
	}
}

