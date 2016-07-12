/**
 * @file src/cstool/main.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <exception>
#include <iostream>
#include <fstream>

#include "../common/xml.hh"
#include "../compile-mat/parser.hh"
#include "cross-section.hh"

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

void usage() {
	std::clog << "Usage: cstool <action> [args]" << std::endl;
	std::clog << "Actions:" << std::endl;
	std::clog << "\tplot cs.xml" << std::endl;
	std::clog << "\tshift cs.xml dx" << std::endl;
	std::clog << "\tmad c1 cs1.xml [c2 cs2.xml]" << std::endl;
	std::clog << "\tmerge cs1.xml cs2.xml x1 x2" << std::endl;
}

int main(int argc, char* argv[]) {
	parser p;
	if(argc < 2) {
		usage();
		return EXIT_SUCCESS;
	}

	std::string action = argv[1];
	if(action == "plot")
	{
		if(argc != 3)
			throw std::runtime_error("plot requires 1 argument");
		std::ifstream ifs1(argv[2]);
		xml::element root1(ifs1);
		std::vector<cstable*> cst_vec;
		for(const xml::element* cst_xe : root1.children("cstable"))
			cst_vec.push_back(cstable::from_xml(*cst_xe));
		generate_tcs_plot(cst_vec, std::cout);
		for(const cstable* cst : cst_vec)
			delete cst;
	}
	else if(action == "shift")
	{
		if(argc != 4)
			throw std::runtime_error("shift requires 2 arguments");
		std::ifstream ifs1(argv[2]);
		double dx = p.eval(argv[3]);
		xml::element root1(ifs1);
		for(const xml::element* cst_xe : root1.children("cstable")) {
			cstable* cst1 = cstable::from_xml(*cst_xe);
			cstable* result = cstable::shift(*cst1, dx);
			std::cout << result->to_xml() << std::endl;
			delete cst1;
			delete result;
		}
	}
	else if(action == "mad")
	{
		if(argc == 4) {
			double c1 = p.eval(argv[2]);
			std::ifstream ifs1(argv[3]);
			xml::element root1(ifs1);
			for(const xml::element* cst_xe : root1.children("cstable")) {
				cstable* cst1 = cstable::from_xml(*cst_xe);
				cstable* result = cstable::mul(c1, *cst1);
				std::cout << result->to_xml() << std::endl;
				delete cst1;
				delete result;
			}
		} else if(argc == 6) {
			double c1 = p.eval(argv[2]);
			std::ifstream ifs1(argv[3]);
			xml::element root1(ifs1);
			double c2 = p.eval(argv[4]);
			std::ifstream ifs2(argv[5]);
			xml::element root2(ifs2);

			if(root1.children("cstable").size() != root2.children("cstable").size())
				throw std::runtime_error("cstable count mismatch");
			uint table_count = root1.children("cstable").size();
			for(uint i=0; i<table_count; i++) {
				cstable* cst1 = cstable::from_xml(*(root1.children("cstable")[i]));
				cstable* cst2 = cstable::from_xml(*(root2.children("cstable")[i]));
				cstable* result = cstable::mad(c1, *cst1, c2, *cst2);
				std::cout << result->to_xml() << std::endl;
				delete cst1;
				delete cst2;
				delete result;
			}
		} else {
			throw std::runtime_error("mad requires 4 or 6 arguments");
		}
	}
	else if(action == "merge")
	{
		if(argc != 6)
			throw std::runtime_error("merge requires 4 arguments");
		std::ifstream ifs1(argv[2]);
		std::ifstream ifs2(argv[3]);
		double e1 = p.eval(argv[4]);
		double e2 = p.eval(argv[5]);
		cstable* tcs1 = cstable::from_xml(*(xml::element(ifs1).children("cstable")[0]));
		cstable* tcs2 = cstable::from_xml(*(xml::element(ifs2).children("cstable")[0]));
		cstable* tcs3 = cstable::merge(*tcs1, *tcs2, e1, e2);
		std::cout << tcs3->to_xml() << std::endl;
		delete tcs1;
		delete tcs2;
		delete tcs3;
	}
	else
	{
		usage();
	}
	return EXIT_SUCCESS;
}
