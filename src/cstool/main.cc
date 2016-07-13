/**
 * @file src/cstool/main.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <exception>
#include <iostream>
#include <fstream>
#include <gsl.h>

#include "../common/xml.hh"
#include "../common/argparse.hh"
#include "../common/command.hh"
#include "../common/parser.hh"

#include "cross-section.hh"

using namespace argparse;
using namespace command;

uint gnuplot_plot_i=0;

void generate_tcs_plot(const std::vector<cstable*> cst_vec, std::ostream& os) {
	os << "set terminal wxt " << gnuplot_plot_i << " persist enhanced" << std::endl;
	os << "set title 'cross section'" << std::endl;
	os << "set xlabel 'K [eV]'" << std::endl;
	os << "set ylabel 'σ [Å^2]'" << std::endl;
	os << "set format x '10^{%T}'" << std::endl;
	os << "set format y '10^{%T}'" << std::endl;
	os << "set logscale xy" << std::endl;
	os << "set grid" << std::endl;
	os << "plot '-' using 1:2 with lines notitle";
	for(uint i=0; i<cst_vec.size()-1; i++)
		os << ", '' using 1:2 with lines notitle";
	os << std::endl;
	for(const cstable* cst : cst_vec) {
		tcstable tcst = cst->to_tcstable();
		for(const std::pair<double,tcs>& tcs_pair : tcst.values) {
			os << (tcs_pair.first/constant::ec);
			os << ' ' << (tcs_pair.second/1e-20);
			os << std::endl;
		}
		os << "end" << std::endl;
	}
	gnuplot_plot_i++;
}

Command cmd_plot("plot",
    "Plot cross-sections in XML file using GnuPlot.",
    [] (gsl::span<std::string> const &argv)
{
    Args args = { 
        positional("<filename>", "XML file for which to plot values.") };
    args.parse(argv);

    std::ifstream ifs1(*args.get<std::string>("<filename>"));

    xml::element root1(ifs1);
    std::vector<cstable*> cst_vec;
    for(const xml::element* cst_xe : root1.children("cstable"))
        cst_vec.push_back(cstable::from_xml(*cst_xe));
    generate_tcs_plot(cst_vec, std::cout);
    for(const cstable* cst : cst_vec)
        delete cst;

    return 0;
});


Command cmd_shift("shift",
    "Shift things.",
    [] (gsl::span<std::string> const &argv)
{ 
    Args args = { 
        positional("<filename>", "XML file for which to shift values."),
        positional("<dx>", "Shift width in units.") };
    args.parse(argv);

    parser p;
    std::ifstream ifs1(*args.get<std::string>("<filename>"));
    double dx = p.eval(*args.get<std::string>("<dx>"));

    xml::element root1(ifs1);
    for(const xml::element* cst_xe : root1.children("cstable")) {
        cstable* cst1 = cstable::from_xml(*cst_xe);
        cstable* result = cstable::shift(*cst1, dx);
        std::cout << result->to_xml() << std::endl;
        delete cst1;
        delete result;
    }

    return 0;
});


Command cmd_mad("mad",
    "Go crazy man.",
    [] (gsl::span<std::string> const &argv)
{ 
    Args args2 = {
        positional("<c1>", "Constant"),
        positional("<filename1>", "XML file") };

    Args args4 = {
        positional("<c1>", "Constant"),
        positional("<filename1>", "XML file"),
        positional("<c2>", "Constant"),
        positional("<filename2>", "XML file") };

    if (argv.size() == 2)
    {
        args2.parse(argv);

        parser p;
        std::ifstream ifs1(*args2.get<std::string>("<filename1>"));
        double c1 = p.eval(*args2.get<std::string>("<c1>"));

        xml::element root1(ifs1);
        for(const xml::element* cst_xe : root1.children("cstable")) 
        {
            cstable* cst1 = cstable::from_xml(*cst_xe);
            cstable* result = cstable::mul(c1, *cst1);
            std::cout << result->to_xml() << std::endl;
            delete cst1;
            delete result;
        }
    }
    else
    {
        args4.parse(argv);

        parser p;
        std::ifstream ifs1(*args4.get<std::string>("<filename1>"));
        double c1 = p.eval(*args4.get<std::string>("<c1>"));
        std::ifstream ifs2(*args4.get<std::string>("<filename2>"));
        double c2 = p.eval(*args4.get<std::string>("<c1>"));

        xml::element root1(ifs1);
        xml::element root2(ifs2);

        if(root1.children("cstable").size() != root2.children("cstable").size())
            throw std::runtime_error("cstable count mismatch");

        uint table_count = root1.children("cstable").size();
        for(uint i=0; i<table_count; i++) 
        {
            cstable* cst1 = cstable::from_xml(*(root1.children("cstable")[i]));
            cstable* cst2 = cstable::from_xml(*(root2.children("cstable")[i]));
            cstable* result = cstable::mad(c1, *cst1, c2, *cst2);
            std::cout << result->to_xml() << std::endl;
            delete cst1;
            delete cst2;
            delete result;
        }
    }

    return 0;
});


Command cmd_merge("merge",
    "Merge cross-section tables.",
    [] (gsl::span<std::string> const &argv)
{
    Args args = {
        positional("<filename1>", "First XML file."),
        positional("<filename2>", "Second XML file."),
        positional("<x1>", "x1"),
        positional("<x2>", "x2") };
    args.parse(argv);

    parser p;
    std::ifstream ifs1(*args.get<std::string>("<filename1>"));
    std::ifstream ifs2(*args.get<std::string>("<filename2>"));
    double e1 = p.eval(*args.get<std::string>("<x1>"));
    double e2 = p.eval(*args.get<std::string>("<x2>"));

    cstable* tcs1 = cstable::from_xml(*(xml::element(ifs1).children("cstable")[0]));
    cstable* tcs2 = cstable::from_xml(*(xml::element(ifs2).children("cstable")[0]));
    cstable* tcs3 = cstable::merge(*tcs1, *tcs2, e1, e2);
    std::cout << tcs3->to_xml() << std::endl;
    delete tcs1;
    delete tcs2;
    delete tcs3;

    return 0;
});

void usage() {
	std::clog << "Usage: cstool <action> [args]" << std::endl;
	std::clog << "Actions:" << std::endl;
	std::clog << "\tplot cs.xml" << std::endl;
	std::clog << "\tshift cs.xml dx" << std::endl;
	std::clog << "\tmad c1 cs1.xml [c2 cs2.xml]" << std::endl;
	std::clog << "\tmerge cs1.xml cs2.xml x1 x2" << std::endl;
}

int main(int argc_, char* argv_[]) {
    std::vector<std::string> argv(argv_ + 1, argv_ + argc_);

    try
    {
        if (argv.size() > 0)
        {
            if (Command::exists(argv[0]))
            {
                Command::at(argv[0])(gsl::as_span(argv).subspan(1));
                return 0;
            } 
        }

        Command::at("list")(argv);
        return 1;
    }

    catch (argparse::Exception const &e)
    {
        std::clog << "Error parsing arguments: " << e.what() << std::endl;
        usage();
        return 1;
    }
}

