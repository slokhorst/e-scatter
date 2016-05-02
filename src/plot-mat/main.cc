/**
 * @file src/plot-mat/main.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <iostream>
#include <fstream>
#include <common/spline.hh>
#include <common/constant.hh>
#include <cdsem/material.hh>

uint gnuplot_plot_i=0;

void generate_elastic_tcs_plot(const material& mat, const double& max_K, std::ostream& os) {
	os << "set terminal wxt " << gnuplot_plot_i << " persist enhanced" << std::endl;
	os << "set title '" << mat.name() << ": elastic cross section'" << std::endl;
	os << "set xlabel 'energy [eV]'" << std::endl;
	os << "set ylabel 'elastic tcs [Å^2]'" << std::endl;
	os << "set format x '10^%T'" << std::endl;
	os << "set format y '10^%T'" << std::endl;
	os << "set logscale xy" << std::endl;
	os << "set grid" << std::endl;
	os << "plot '-' using 1:2 with lines notitle" << std::endl;
	os << "# kinetic energy [eV] elastic tcs [Å^2]" << std::endl;
	for(double K = mat.fermi(); K < max_K; K += constant::ec) {
		os << (K/constant::ec);
		os << ' ' << (mat.elastic_tcs(K)/1e-20);
		os << std::endl;
	}
	os << "end" << std::endl;
	gnuplot_plot_i++;
}

void generate_inelastic_tcs_plot(const material& mat, const double& max_K, std::ostream& os) {
	os << "set terminal wxt " << gnuplot_plot_i << " persist enhanced" << std::endl;
	os << "set title '" << mat.name() << ": inelastic cross section'" << std::endl;
	os << "set xlabel 'energy [eV]'" << std::endl;
	os << "set ylabel 'inelastic tcs [Å^2]'" << std::endl;
	os << "set format x '10^%T'" << std::endl;
	os << "set format y '10^%T'" << std::endl;
	os << "set logscale xy" << std::endl;
	os << "set grid" << std::endl;
	os << "plot '-' using 1:2 with lines notitle" << std::endl;
	os << "# kinetic energy [eV] inelastic tcs [Å^2]" << std::endl;
	for(double K = mat.fermi(); K < max_K; K += constant::ec) {
		os << (K/constant::ec);
		os << ' ' << (mat.inelastic_tcs(K)/1e-20);
		os << std::endl;
	}
	os << "end" << std::endl;
	gnuplot_plot_i++;
}

void generate_elastic_mfp_plot(const material& mat, const double& max_K, std::ostream& os) {
	os << "set terminal wxt " << gnuplot_plot_i << " persist enhanced" << std::endl;
	os << "set title '" << mat.name() << ": elastic mean free path'" << std::endl;
	os << "set xlabel 'energy [eV]'" << std::endl;
	os << "set ylabel 'elastic mfp [nm]'" << std::endl;
	os << "set format x '10^%T'" << std::endl;
	os << "set format y '10^%T'" << std::endl;
	os << "set logscale xy" << std::endl;
	os << "set grid" << std::endl;
	os << "plot '-' using 1:2 with lines notitle" << std::endl;
	os << "# kinetic energy [eV] elastic mfp [nm]" << std::endl;
	for(double K = mat.fermi(); K < max_K; K += constant::ec) {
		os << (K/constant::ec);
		os << ' ' << (1/(mat.density()*mat.elastic_tcs(K))/1e-9);
		os << std::endl;
	}
	os << "end" << std::endl;
	gnuplot_plot_i++;
}

void generate_inelastic_mfp_plot(const material& mat, const double& max_K, std::ostream& os) {
	os << "set terminal wxt " << gnuplot_plot_i << " persist enhanced" << std::endl;
	os << "set title '" << mat.name() << ": inelastic mean free path'" << std::endl;
	os << "set xlabel 'energy [eV]'" << std::endl;
	os << "set ylabel 'inelastic mfp [nm]'" << std::endl;
	os << "set format x '10^%T'" << std::endl;
	os << "set format y '10^%T'" << std::endl;
	os << "set logscale xy" << std::endl;
	os << "set grid" << std::endl;
	os << "plot '-' using 1:2 with lines notitle" << std::endl;
	os << "# kinetic energy [eV] inelastic mfp [nm]" << std::endl;
	for(double K = mat.fermi(); K < max_K; K += constant::ec) {
		os << (K/constant::ec);
		os << ' ' << (1/(mat.density()*mat.inelastic_tcs(K))/1e-9);
		os << std::endl;
	}
	os << "end" << std::endl;
	gnuplot_plot_i++;
}

void generate_elastic_scatterangle_distribution(const material& mat, const double& max_K, std::ostream& os) {
	os << "set terminal wxt " << gnuplot_plot_i << " persist enhanced" << std::endl;
	os << "set title '" << mat.name() << ": mean elastic scatter angle'" << std::endl;
	os << "set xlabel 'energy [eV]'" << std::endl;
	os << "set ylabel 'mean angle'" << std::endl;
	os << "set yrange [0:pi]" << std::endl;
	os << "set format x '10^%T'" << std::endl;
	os << "unset format y" << std::endl;
	os << "set ytics (0, 'π/2' 0.5*pi, 'π' pi)" << std::endl;
	os << "set logscale x" << std::endl;
	os << "unset logscale y" << std::endl;
	os << "set grid" << std::endl;
	os << "plot '-' using 1:2:3:4 with yerrorbars notitle" << std::endl;
	os << "# kinetic energy [eV] mean angle [rad]" << std::endl;
	for(const auto& energy_p_dcs_map : mat._elastic_dcs){
		const double& K = std::exp(energy_p_dcs_map.first);
		std::map<double,double> theta_p_map;
		for(const auto& p_dcs_pair : energy_p_dcs_map.second)
			theta_p_map[p_dcs_pair.second] = p_dcs_pair.first;
		spline cdf = spline::linear(theta_p_map);
		const double mean = constant::pi - cdf.integrate(0)(constant::pi);
		const double lower = mat.elastic_dcs(K,0.05);
		const double upper = mat.elastic_dcs(K,0.95);

		os << (K/constant::ec);
		os << ' ' << (mean);
		os << ' ' << (lower);
		os << ' ' << (upper);
		os << std::endl;
	}
	os << "end" << std::endl;
	gnuplot_plot_i++;
}

void generate_inelastic_energyloss_distribution(const material& mat, const double& max_K, std::ostream& os) {
	os << "set terminal wxt " << gnuplot_plot_i << " persist enhanced" << std::endl;
	os << "set title '" << mat.name() << ": mean inelastic energy loss'" << std::endl;
	os << "set xlabel 'energy [eV]'" << std::endl;
	os << "set ylabel 'energy loss [eV]'" << std::endl;
	os << "set format x '10^%T'" << std::endl;
	os << "set format y '10^%T'" << std::endl;
	os << "set logscale xy" << std::endl;
	os << "set grid" << std::endl;
	os << "plot '-' using 1:2:3:4 with yerrorbars notitle" << std::endl;
	os << "# kinetic energy [eV] mean energy loss [eV]" << std::endl;
	for(const auto& energy_p_dcs_map : mat._inelastic_dcs){
		const double& K = std::exp(energy_p_dcs_map.first);
		std::map<double,double> omega_p_map;
		for(const auto& p_dcs_pair : energy_p_dcs_map.second)
			omega_p_map[p_dcs_pair.second] = p_dcs_pair.first;
		spline cdf = spline::linear(omega_p_map);
		double omega2 = K-mat.fermi();
		const double mean = omega2 - cdf.integrate(0)(omega2);
		const double lower = mat.inelastic_dcs(K,0.05);
		const double upper = mat.inelastic_dcs(K,0.95);

		os << (K/constant::ec);
		os << ' ' << (mean/constant::ec);
		os << ' ' << (lower/constant::ec);
		os << ' ' << (upper/constant::ec);
		os << std::endl;
	}
	os << "end" << std::endl;
	gnuplot_plot_i++;
}

void generate_ionization_plot(const material& mat, const double& max_K, std::ostream& os) {
	os << "set terminal wxt " << gnuplot_plot_i << " persist enhanced" << std::endl;
	os << "set title '" << mat.name() << ": ionization mean free path'" << std::endl;
	os << "set xlabel 'energy [eV]'" << std::endl;
	os << "set ylabel 'ionization mfp [nm]'" << std::endl;
	os << "set format x '10^%T'" << std::endl;
	os << "unset format y" << std::endl;
	os << "set logscale x" << std::endl;
	os << "unset logscale y" << std::endl;
	os << "set grid" << std::endl;
	os << "plot '-' using 1:2:3:4 with yerrorbars notitle" << std::endl;
	os << "# kinetic energy [eV] ionization mfp [nm]" << std::endl;
	for(double K = mat.fermi(); K < max_K; K += 0.1*constant::ec) {
		os << (K/constant::ec);
		os << ' ' << (mat.ionization_energy(K,0.5)/constant::ec);
		os << ' ' << (mat.ionization_energy(K,0.05)/constant::ec);
		os << ' ' << (mat.ionization_energy(K,0.95)/constant::ec);
		os << std::endl;
	}
	os << "end" << std::endl;
	gnuplot_plot_i++;
}

int main(int argc, char* argv[]) {
	std::string input_file = argv[1];
	std::ifstream ifs(input_file);
	archive::istream is(ifs);

	material mat("tmp", 0, 0, 0);
	is >> mat;

	generate_elastic_mfp_plot(mat, 1e4*constant::ec, std::cout);
	//generate_elastic_scatterangle_distribution(mat, 1e4*constant::ec, std::cout);

	generate_inelastic_mfp_plot(mat, 1e4*constant::ec, std::cout);
	//generate_inelastic_energyloss_distribution(mat, 1e4*constant::ec, std::cout);

	generate_ionization_plot(mat, 1e4*constant::ec, std::cout);
}
