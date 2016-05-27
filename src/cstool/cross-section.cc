/**
 * @file src/cstool/cross-section.cc
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <cstool/cross-section.hh>
#include <set>
#include <common/interpolate.hh>
#include <common/spline.hh>
#include <cstool/cross-section.hh>
#include <compile-mat/parser.hh>

//helper functions

std::pair<double,tcs> tcspair_from_xml(const std::string& prop, const xml::element& parent) {
	parser p;
	const double value = p.eval(parent.attr(prop));
	const double tcs = p.eval(parent.attr("cs"));
	return std::make_pair(value, tcs);
}
std::string tcspair_to_xml(const std::string& prop, const std::pair<double,tcs>& tcs_pair) {
	std::stringstream xmlss;
	xmlss << "\t<cross-section " << prop << "=\"" << tcs_pair.first << "\" cs=\"" << tcs_pair.second << "\" />";
	return xmlss.str();
}

std::pair<double,dcs> dcspair_from_xml(const std::string& prop, const std::string& diff_prop, const xml::element& parent) {
	parser p;
	const double value = p.eval(parent.attr(prop));
	dcs dcs;
	for(const xml::element* dcs_p : parent.children("insert")) {
		double diff_value = p.eval(dcs_p->attr(diff_prop));
		double dcs_value = p.eval(dcs_p->attr("dcs"));
		dcs[diff_value] = dcs_value;
	}
	return std::make_pair(value, dcs);
}
std::string dcspair_to_xml(const std::string& prop, const std::string& diff_prop, const std::pair<double,dcs>& dcs_pair) {
	std::stringstream xmlss;
	xmlss << "\t<cross-section " << prop << "=\"" << dcs_pair.first << "\">" << std::endl;
	for(const std::pair<double,double>& cs : dcs_pair.second)
		xmlss << "\t\t<insert " << diff_prop << "=\"" << cs.first << "\" dcs=\"" << cs.second << "\" />" << std::endl;
	xmlss << "\t</cross-section>";
	return xmlss.str();
}

tcs dcs_to_tcs(const std::string& prop, const std::string& diff_prop, const std::function<double(double)>& differential, const std::pair<double,double>& int_limits, const dcs& dcs) {
	double tcs;
	std::map<double,double> dcs_int_map;
	for(const std::pair<double,double>& d : dcs) {
		const double& x = d.first;
		const double& y = d.second;
		if((x >= int_limits.first) && (x <= int_limits.second) && (y >= 0))
			dcs_int_map[x] = y*differential(x);
	}
	if(dcs_int_map.empty()) {
		tcs = 0;
	} else {
		dcs_int_map[int_limits.first] = 0;
		dcs_int_map[int_limits.second] = 0;
		tcs = spline::linear(dcs_int_map).integrate(0)(int_limits.second);
	}

	return tcs;
}

//cstable

std::string cstable::type() const {
	return m_type;
}
std::string cstable::attribute(const std::string& attr) const {
	return m_attr_map.at(attr);
}

cstable* cstable::from_xml(const xml::element& parent) {
	std::string type = parent.attr("type");
	if (type == "elastic" || type == "inelastic")
		return dcstable::from_xml(parent);
	else
		return tcstable::from_xml(parent);
}
cstable* cstable::mul(double c_a, const cstable& a) {
	if(typeid(a) == typeid(tcstable)) {
		const tcstable tcst = dynamic_cast<const tcstable&>(a);
		return tcstable::mul(c_a, tcst);
	} else if(typeid(a) == typeid(dcstable)) {
		const dcstable dcst = dynamic_cast<const dcstable&>(a);
		return dcstable::mul(c_a, dcst);
	} else {
		throw std::runtime_error("unknown cstable");
	}
}
cstable* cstable::mad(double c_a, const cstable& a, double c_b, const cstable& b) {
	if(typeid(a) == typeid(tcstable) && typeid(b) == typeid(tcstable)) {
		const tcstable tcst1 = dynamic_cast<const tcstable&>(a);
		const tcstable tcst2 = dynamic_cast<const tcstable&>(b);
		return tcstable::mad(c_a, tcst1, c_b, tcst2);
	} else if(typeid(a) == typeid(dcstable) && typeid(b) == typeid(dcstable)) {
		const dcstable dcst1 = dynamic_cast<const dcstable&>(a);
		const dcstable dcst2 = dynamic_cast<const dcstable&>(b);
		return dcstable::mad(c_a, dcst1, c_b, dcst2);
	} else {
		throw std::runtime_error("unknown cstable");
	}
}
cstable* cstable::merge(const cstable& a, const cstable& b, double x1, double x2) {
	if(typeid(a) == typeid(tcstable) && typeid(b) == typeid(tcstable)) {
		const tcstable tcst1 = dynamic_cast<const tcstable&>(a);
		const tcstable tcst2 = dynamic_cast<const tcstable&>(b);
		return tcstable::merge(tcst1, tcst2, x1, x2);
	} else if(typeid(a) == typeid(dcstable) && typeid(b) == typeid(dcstable)) {
		const dcstable dcst1 = dynamic_cast<const dcstable&>(a);
		const dcstable dcst2 = dynamic_cast<const dcstable&>(b);
		return dcstable::merge(dcst1, dcst2, x1, x2);
	} else {
		throw std::runtime_error("unknown cstable");
	}
}

//tcstable

std::string tcstable::prop() const {
	return m_type_prop_map.at(m_type);
}
tcstable tcstable::to_tcstable() const {
	tcstable tcst(*this);
	return tcst;
}
std::string tcstable::to_xml() const {
	std::stringstream xmlss;
	xmlss << "<cstable type=\""+type()+"\"";
	for(const std::pair<std::string,std::string> pair : m_attr_map)
		xmlss << " " << pair.first << "=\"" << pair.second << "\"";
	xmlss << ">" << std::endl;
	for(const std::pair<double,tcs>& cs : values){
		xmlss << tcspair_to_xml(prop(), cs) << std::endl;
	}
	xmlss << "</cstable>" << std::endl;
	return xmlss.str();
}
tcstable* tcstable::from_xml(const xml::element& parent) {
	parser p;
	tcstable* tcst = new tcstable(parent.attr("type"));
	//TODO: check if required attributes are present
	//tcst->m_attr_map = parent.attributes();
	tcst->m_attr_map.erase("type");
	for(const xml::element* cross_section_p : parent.children("cross-section"))
		tcst->values.insert(tcspair_from_xml(tcst->prop(), *cross_section_p));
	if(tcst->m_attr_map.count("occupancy") > 0){
		double occupancy = p.eval(tcst->m_attr_map.at("occupancy"));
		tcst->m_attr_map.erase("occupancy");
		tcstable* tcst2 = tcstable::mul(occupancy,*tcst);
		delete tcst;
		tcst = tcst2;
	}
	return tcst;
}
tcstable* tcstable::mul(double c_a, const tcstable& a) {
	tcstable* res = new tcstable(a);
	for(auto& pair : res->values)
		pair.second *= c_a;
	return res;
}
tcstable* tcstable::mad(double c_a, const tcstable& a, double c_b, const tcstable& b) {
	if(a.type() != b.type())
		throw std::runtime_error("types must match");
	tcstable* result = new tcstable(a.type());
	std::set<double> x_union;
	for(const auto pair : a.values)
		x_union.insert(pair.first);
	for(const auto pair : b.values)
		x_union.insert(pair.first);
	for(const double x : x_union)
		result->values[x] = c_a*interpolate(a.values,x) + c_b*interpolate(b.values,x);
	return result;
}
tcstable* tcstable::merge(const tcstable& a, const tcstable& b, double x1, double x2) {
	if(a.type() != b.type())
		throw std::runtime_error("types must match");
	tcstable* result = new tcstable(a.type());
	std::set<double> x_merge_points;
	for(const auto& a1 : a.values)
		if(a1.first <= x1) {
			result->values[a1.first] = a1.second;
		} else if(a1.first <= x2) {
			x_merge_points.insert(a1.first);
		}
	for(const auto& b1 : b.values)
		if(b1.first >= x2) {
			result->values[b1.first] = b1.second;
		} else if(b1.first >= x1) {
			x_merge_points.insert(b1.first);
		}
	x_merge_points.insert(x1);
	x_merge_points.insert(x2);
	for(const double& x : x_merge_points) {
		const double t = (std::log(x)-std::log(x1))/(std::log(x2)-std::log(x1));
		result->values[x] = (1.0-t)*interpolate(a.values,x) + t*interpolate(b.values,x);
	}
	return result;
}


//dcstable

std::string dcstable::prop() const {
	return m_type_prop_map.at(m_type).first;
}
std::string dcstable::diff_prop() const {
	return m_type_prop_map.at(m_type).second;
}
std::function<double(double)> dcstable::differential() const {
	return m_type_differential_map.at(m_type);
}
tcstable dcstable::to_tcstable() const {
	tcstable tcst(type());
	for(const std::pair<double,dcs>& dcs_pair : values) {
		std::pair<double,double> int_limits;
		if(m_type == "elastic")
			int_limits = std::make_pair(0,constant::pi);
		else
			int_limits = std::make_pair(0,dcs_pair.second.crbegin()->first);
		tcst.values[dcs_pair.first] = dcs_to_tcs(prop(), diff_prop(), differential(), int_limits, dcs_pair.second);
	}
	return tcst;
}
std::string dcstable::to_xml() const {
	std::stringstream xmlss;
	xmlss << "<cstable type=\""+type()+"\"";
	for(const std::pair<std::string,std::string> pair : m_attr_map)
		xmlss << " " << pair.first << "=\"" << pair.second << "\"";
	xmlss << ">" << std::endl;
	for(const std::pair<double,dcs>& dcs_pair : values)
		xmlss << dcspair_to_xml(prop(), diff_prop(), dcs_pair) << std::endl;
	xmlss << "</cstable>" << std::endl;
	return xmlss.str();
}
dcstable* dcstable::from_xml(const xml::element& parent) {
	parser p;
	dcstable* dcst = new dcstable(parent.attr("type"));
	//TODO: check if required attributes are present
	//dcst->m_attr_map = parent.attributes();
	dcst->m_attr_map.erase("type");
	for(const xml::element* cross_section_p : parent.children("cross-section"))
		dcst->values.insert(dcspair_from_xml(dcst->prop(), dcst->diff_prop(), *cross_section_p));
	if(dcst->m_attr_map.count("occupancy") > 0){
		double occupancy = p.eval(dcst->m_attr_map.at("occupancy"));
		dcst->m_attr_map.erase("occupancy");
		dcstable* dcst2 = dcstable::mul(occupancy,*dcst);
		delete dcst;
		dcst = dcst2;
	}
	return dcst;
}
dcstable* dcstable::mul(double c_a, const dcstable& a) {
	dcstable* res = new dcstable(a);
	for(auto& pair : res->values)
		for(auto& pair2 : pair.second)
			pair2.second *= c_a;
	return res;
}
dcstable* dcstable::mad(double c_a, const dcstable& a, double c_b, const dcstable& b) {
	if(a.type() != b.type())
		throw std::runtime_error("types must match");
	dcstable* result = new dcstable(a.type());
	std::set<double> x_union;
	for(const auto pair : a.values)
		x_union.insert(pair.first);
	for(const auto pair : b.values)
		x_union.insert(pair.first);
	for(const double x : x_union) {
		std::set<double> y_union;
		if(a.values.count(x))
			for(const auto pair : a.values.at(x))
				y_union.insert(pair.first);
		if(b.values.count(x))
			for(const auto pair : b.values.at(x))
				y_union.insert(pair.first);
		for(const double y : y_union)
			result->values[x][y] = c_a*interpolate(a.values,x,y) + c_b*interpolate(b.values,x,y);
	}
	return result;
}
dcstable* dcstable::merge(const dcstable& a, const dcstable& b, double x1, double x2) {
	if(a.type() != b.type())
		throw std::runtime_error("types must match");
	dcstable* result = new dcstable(a.type());
	std::set<double> x_merge_points, y_merge_points;
	for(const auto& a1 : a.values)
		if(a1.first <= x1) {
			result->values[a1.first] = a1.second;
		} else if(a1.first <= x2) {
			x_merge_points.insert(a1.first);
			for(const auto& a2 : a1.second)
				y_merge_points.insert(a2.first);
		}
	for(const auto& b1 : b.values)
		if(b1.first >= x2) {
			result->values[b1.first] = b1.second;
		} else if(b1.first >= x1) {
			x_merge_points.insert(b1.first);
			for(const auto& b2 : b1.second)
				y_merge_points.insert(b2.first);
		}
	x_merge_points.insert(x1);
	x_merge_points.insert(x2);
	for(const double& x : x_merge_points) {
		const double t = (std::log(x)-std::log(x1))/(std::log(x2)-std::log(x1));
		for(const double& y : y_merge_points)
			result->values[x][y] = (1.0-t)*interpolate(a.values,x,y) + t*interpolate(b.values,x,y);
	}
	return result;
}
