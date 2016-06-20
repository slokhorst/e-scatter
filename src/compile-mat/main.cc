/**
 * @file src/compile-mat/main.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <common/constant.hh>
#include <common/archive.hh>
#include <common/optional.hh>
#include <common/xml.hh>
#include <compile-mat/parser.hh>
#include <cdsem/material.hh>
#include <cstool/cross-section.hh>

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
	po::options_description opts("allowed options");
	opts.add_options()
		("help,h", "produce help message")
		("output", po::value<std::string>()->required(), "output filename for binary material (.mat)")
		("name", po::value<std::string>()->required(), "material name")
		("elastic", po::value<std::string>()->required(), "input filename for elastic cross-sections xml")
		("inelastic", po::value<std::string>()->required(), "input filename for inelastic cross-sections xml")
		("ionization", po::value<std::string>()->required(), "input filename for ionization cross-sections xml")
		("number-density", po::value<std::string>()->required(), "number density [#/m^3]")
		("fermi-energy", po::value<std::string>()->required(), "Fermi energy [J]")
		("work-function", po::value<std::string>()->required(), "work function [J]")
		("band-gap", po::value<std::string>(), "band gap [J]")
		("phonon-loss", po::value<std::string>(), "phonon loss [J]");
	po::positional_options_description pos_opts;
	pos_opts.add("output", 1);
	po::variables_map var_map;
	try {
		po::store(po::command_line_parser(argc, argv).options(opts).positional(pos_opts).run(), var_map);
		if(var_map.count("help") > 0) {
			std::clog << opts;
			return EXIT_SUCCESS;
		}
		po::notify(var_map);
	} catch(const std::exception& e) {
		std::clog << "error parsing program options: " << e.what() << std::endl;
		std::clog << opts;
		return EXIT_FAILURE;
	}

	std::string output_fn = var_map["output"].as<std::string>();

	try {
		parser p;
		std::string name = var_map["name"].as<std::string>();
		double number_density = p.eval(var_map["number-density"].as<std::string>());
		double fermi = p.eval(var_map["fermi-energy"].as<std::string>());
		double work_func = p.eval(var_map["work-function"].as<std::string>());
		optional<double> bandgap;
		if(var_map.count("band-gap") > 0)
			bandgap = p.eval(var_map["band-gap"].as<std::string>());
		double phononloss = p.eval(var_map["phonon-loss"].as<std::string>());

		material mat(name, fermi, (fermi+work_func), phononloss, number_density);
		if(bandgap.is_defined())
			mat = material(name, fermi, (fermi+work_func), phononloss, bandgap(), number_density);

		std::ifstream elastic_ifs(var_map["elastic"].as<std::string>());
		xml::element elastic_root(elastic_ifs);
		if(elastic_root.children("cstable").size()!=1)
			throw std::runtime_error("need exactly one elastic dcstable");
		dcstable* el_dcst = dcstable::from_xml(*(elastic_root.children("cstable")[0]));
		for(auto energy_data : el_dcst->values)
			mat.set_elastic_data(energy_data.first, energy_data.second);
		delete el_dcst;

		std::ifstream inelastic_ifs(var_map["inelastic"].as<std::string>());
		xml::element inelastic_root(inelastic_ifs);
		if(inelastic_root.children("cstable").size()!=1)
			throw std::runtime_error("need exactly one inelastic dcstable");
		dcstable* inel_dcst = dcstable::from_xml(*(inelastic_root.children("cstable")[0]));
		for(auto energy_data : inel_dcst->values)
			mat.set_inelastic_data(energy_data.first, energy_data.second);
		delete inel_dcst;

		std::ifstream ionization_ifs(var_map["ionization"].as<std::string>());
		xml::element ionization_root(ionization_ifs);
		if(ionization_root.children("cstable").size()<1)
			std::clog << "WARNING: no ionization cross-sections";
		for(const xml::element* ion_xml_p : ionization_root.children("cstable")) {
			tcstable* tcst = tcstable::from_xml(*ion_xml_p);
			double B = p.eval(tcst->attribute("binding-energy"));
			std::map<double,tcs> values_plus_B;
			for(const std::pair<double,tcs>& tcs_pair : tcst->values)
				values_plus_B[tcs_pair.first+B] = tcs_pair.second;
			mat.set_ionization_data(B,values_plus_B);
			delete tcst;
		}

		std::clog  << "material name = `" << mat.name() << "`"
			<< ", number-density = " << mat.density() << " [#/m^3]"
			<< ", fermi-energy = " << (mat.fermi()/constant::ec) << " [eV]"
			<< ", barrier = " << (mat.barrier()/constant::ec) << " [eV]"
			<< ", phonon-loss = " << (mat.phononloss()/constant::ec) << " [eV]";
		if(mat.bandgap().is_defined())
			std::clog << ", band-gap = " << (mat.bandgap()()/constant::ec) << " [eV]";
		std::clog << std::endl;

		std::clog << "serializing material to `" << output_fn << "`" << std::endl;
		std::ofstream ofs(output_fn);
		if(!ofs)
			throw std::ios_base::failure("ofstream open failure");
		archive::ostream os(ofs);
		os << mat;
	} catch(const std::exception& e) {
		std::clog << "failed to serialize to `" << output_fn << "`"
			<< "\n\terror: " << e.what() << std::endl;
	}
	std::clog << "serialization completed" << std::endl;
	return EXIT_SUCCESS;
}
