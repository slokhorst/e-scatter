/**
 * @file src/compile-mat/main.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cpl/mathexpr.h>
#include <cpl/text.h>
#include <common/constant.hh>
#include <cpl/xml.h>
#include <cdsem/material.hh>
#include <cstool/cross-section.hh>
#include <common/archive.hh>
#include <common/optional.hh>

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
		("band-gap", po::value<std::string>(), "band gap [J]");
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
		std::clog << cpl::text::cat("error parsing program options: ",e.what()) << std::endl;
		std::clog << opts;
		return EXIT_FAILURE;
	}

	std::string output_fn = var_map["output"].as<std::string>();

	try {
		std::string name = var_map["name"].as<std::string>();
		double number_density = cpl::mathexpr(var_map["number-density"].as<std::string>()).evaluate();
		double fermi = cpl::mathexpr(var_map["fermi-energy"].as<std::string>()).evaluate();
		double work_func = cpl::mathexpr(var_map["work-function"].as<std::string>()).evaluate();
		optional<double> bandgap;
		if(var_map.count("band-gap") > 0)
			bandgap = cpl::mathexpr(var_map["band-gap"].as<std::string>()).evaluate();

		material mat(name, fermi, (fermi+work_func), number_density);
		if(bandgap.is_defined())
			mat = material(name, fermi, (fermi+work_func), bandgap(), number_density);

		cpl::xml::element elastic_root(var_map["elastic"].as<std::string>());
		if(elastic_root.children("cstable").size()!=1)
			throw std::runtime_error("need exactly one elastic dcstable");
		dcstable* el_dcst = dcstable::from_xml(*(elastic_root.children("cstable")[0]));
		for(auto energy_data : el_dcst->values)
			mat.set_elastic_data(energy_data.first, energy_data.second);
		delete el_dcst;

		cpl::xml::element inelastic_root(var_map["inelastic"].as<std::string>());
		if(inelastic_root.children("cstable").size()!=1)
			throw std::runtime_error("need exactly one inelastic dcstable");
		dcstable* inel_dcst = dcstable::from_xml(*(inelastic_root.children("cstable")[0]));
		for(auto energy_data : inel_dcst->values)
			mat.set_inelastic_data(energy_data.first, energy_data.second);
		delete inel_dcst;

		cpl::xml::element ionization_root(var_map["ionization"].as<std::string>());
		if(ionization_root.children("cstable").size()<1)
			std::clog << "WARNING: no ionization cross-sections";
		for(const cpl::xml::element* ion_xml_p : ionization_root.children("cstable")) {
			tcstable* tcst = tcstable::from_xml(*ion_xml_p);
			double B = cpl::mathexpr(tcst->attribute("binding-energy")).evaluate();
			std::map<double,tcs> values_plus_B;
			for(const std::pair<double,tcs>& tcs_pair : tcst->values)
				values_plus_B[tcs_pair.first+B] = tcs_pair.second;
			mat.set_ionization_data(B,values_plus_B);
			delete tcst;
		}

		std::string log_msg = cpl::text::cat("material name = `", mat.name(), "`");
		log_msg += cpl::text::cat(", number-density = ", mat.density(), " [#/m^3]");
		log_msg += cpl::text::cat(", fermi-energy = ", mat.fermi()/constant::ec, " [eV]");
		log_msg += cpl::text::cat(", barrier = ", mat.barrier()/constant::ec, " [eV]");
		if(mat.bandgap().is_defined())
			log_msg += cpl::text::cat(", band-gap = ", mat.bandgap()()/constant::ec, " [eV]");
		std::clog << log_msg << std::endl;

		std::cerr << cpl::text::cat("serializing material to `", output_fn, "`") << std::endl;
		std::ofstream ofs(output_fn);
		if(!ofs)
			throw std::ios_base::failure("ofstream open failure");
		archive::ostream os(ofs);
		os << mat;
	} catch(const std::exception& e) {
		std::cerr << cpl::text::cat("failed to serialize to `", output_fn,"`",
			"\n\terror: ",e.what()) << std::endl;
	}
	std::clog << "serialization completed" << std::endl;
	return EXIT_SUCCESS;
}
