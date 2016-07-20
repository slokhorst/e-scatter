/**
 * @file src/cstool/compile-mat.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 * @author Johan Hidding <J.Hidding@esciencecenter.nl>
 */

#include <iostream>
#include <fstream>
#include <experimental/optional>

#include "../common/constant.hh"
#include "../common/archive.hh"
//#include "../common/optional.hh"
#include "../common/xml.hh"
#include "../common/parser.hh"

#include "../common/material.hh"

#include "../common/argparse.hh"
#include "../common/command.hh"

#include "cross-section.hh"

using namespace argparse;
using command::Command;

Command cmd_compile_mat("compile-mat",
    "Compile a material to the binary input files that are used by CdSEM & Co.",
    [] (gsl::span<std::string> const &argv)
{
    Args args = {
        flag("--help", "produce help message"),
        positional("<output>",     "output file for binary material (.mat)"),
        option("--name",           "material name"),
        option("--elastic",        "input XML file for elastic cross-sections"),
        option("--inelastic",      "input XML file for inelastic cross-sections"),
        option("--ionization",     "input XML file for ionization cross-sections"),
        option("--number-density", "number density [#/m^3]"),
        option("--fermi-energy",   "Fermi energy [J]"),
        option("--work-function",  "work function [J]"),
        option("--band-gap",       "band gap [J]", "", true)
    };
    try
    {
        args.parse(argv);
    }
    catch (argparse::Exception const &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << LongString(
            "Compile a material to the binary input files that are used by CdSEM & Co.",
            80, [] () { return fancy::compose(""); }) << "\n\n";
        std::cerr << args.usage();
        return EXIT_FAILURE;
    }

    if (*args.get<bool>("--help"))
    {
        std::cerr << LongString(
            "Help requested. Compile a material to the binary input files that are used by CdSEM & Co.",
            80, [] () { return fancy::compose(""); }) << "\n\n";
        std::cerr << args.usage();
        return EXIT_FAILURE;
    }

    parser p;
	std::string output_fn = *args.get<std::string>("<output>");
    std::string name = *args.get<std::string>("--name");
    double number_density = p.eval(*args.get<std::string>("--number-density"));
    double fermi = p.eval(*args.get<std::string>("--fermi-energy"));
    double work_func = p.eval(*args.get<std::string>("--work-function"));
    std::experimental::optional<double> bandgap = args.get<double>("--band-gap");

    material mat(name, fermi, (fermi+work_func), number_density);
    if (bandgap)
        mat = material(name, fermi, (fermi+work_func), *bandgap, number_density);

    std::ifstream elastic_ifs(*args.get<std::string>("--elastic"));
    xml::element elastic_root(elastic_ifs);
    if(elastic_root.children("cstable").size()!=1)
        throw std::runtime_error("need exactly one elastic dcstable");
    dcstable* el_dcst = dcstable::from_xml(*(elastic_root.children("cstable")[0]));
    for(auto energy_data : el_dcst->values)
        mat.set_elastic_data(energy_data.first, energy_data.second);
    delete el_dcst;

    std::ifstream inelastic_ifs(*args.get<std::string>("--inelastic"));
    xml::element inelastic_root(inelastic_ifs);
    if(inelastic_root.children("cstable").size()!=1)
        throw std::runtime_error("need exactly one inelastic dcstable");
    dcstable* inel_dcst = dcstable::from_xml(*(inelastic_root.children("cstable")[0]));
    for(auto energy_data : inel_dcst->values)
        mat.set_inelastic_data(energy_data.first, energy_data.second);
    delete inel_dcst;

    std::ifstream ionization_ifs(*args.get<std::string>("--ionization"));
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

    std::clog << "material name = `" << mat.name() << "`"
              << ", number-density = " << mat.density() << " [#/m^3]"
              << ", fermi-energy = " << (mat.fermi()/constant::ec) << " [eV]"
              << ", barrier = " << (mat.barrier()/constant::ec) << " [eV]";

    if(mat.bandgap().is_defined())
        std::clog << ", band-gap = " << (mat.bandgap()()/constant::ec) << " [eV]";

    std::clog << std::endl;

    std::clog << "serializing material to `" << output_fn << "`" << std::endl;

    std::ofstream ofs(output_fn);
    if(!ofs)
        throw std::ios_base::failure("ofstream open failure");

    archive::ostream os(ofs);
    os << mat;

	std::clog << "serialization completed" << std::endl;
	return EXIT_SUCCESS;
});
