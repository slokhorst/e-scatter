/**
 * @file src/cstool/cross-section.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__COMMON__CROSS_SECTION__HEADER_INCLUDED
#define eSCATTER__COMMON__CROSS_SECTION__HEADER_INCLUDED

#include <functional>
#include <cmath>
#include <common/constant.hh>
#include <common/xml.hh>

using tcs = double;
using dcs = std::map<double,double>;

class cstable;
class tcstable;
class dcstable;

class cstable {
public:
	cstable(const std::string& t) : m_type(t) {}
	virtual ~cstable() = default;
	virtual std::string type() const;
	virtual std::string attribute(const std::string& attr) const;
	virtual std::string prop() const = 0;
	virtual tcstable to_tcstable() const = 0;
	virtual std::string to_xml() const = 0;
	static cstable* from_xml(const xml::element& parent);
	static cstable* mul(double c_a, const cstable& a);
	static cstable* mad(double c_a, const cstable& a, double c_b, const cstable& b);
	static cstable* merge(const cstable& a, const cstable& b, double x1, double x2);
protected:
	std::string m_type;
	std::map<std::string,std::string> m_attr_map;
	std::map<std::string,std::vector<std::string>> m_type_required_attr_map = {
		{"elastic", {}},
		{"inelastic", {}},
		{"ionization", {"binding-energy"}}
	};
};

class tcstable : public cstable {
public:
	tcstable(const std::string& t) : cstable(t) {}
	~tcstable() = default;
	std::string prop() const;
	tcstable to_tcstable() const;
	std::string to_xml() const;
	static tcstable* from_xml(const xml::element& parent);
	static tcstable* mul(double c_a, const tcstable& a);
	static tcstable* mad(double c_a, const tcstable& a, double c_b, const tcstable& b);
	static tcstable* merge(const tcstable& a, const tcstable& b, double x1, double x2);
	std::map<double,tcs> values;
private:
	std::map<std::string,std::string> m_type_prop_map = {
		{"elastic", "energy"},
		{"inelastic", "energy"},
		{"ionization", "energy"}
	};
};

class dcstable : public cstable {
public:
	dcstable(const std::string& t) : cstable(t) {}
	~dcstable() = default;
	std::string prop() const;
	std::string diff_prop() const;
	std::function<double(double)> differential() const;
	tcstable to_tcstable() const;
	std::string to_xml() const;
	static dcstable* from_xml(const xml::element& parent);
	static dcstable* mul(double c_a, const dcstable& a);
	static dcstable* mad(double c_a, const dcstable& a, double c_b, const dcstable& b);
	static dcstable* merge(const dcstable& a, const dcstable& b, double x1, double x2);
	std::map<double,dcs> values;
private:
	std::map<std::string,std::pair<std::string,std::string>> m_type_prop_map = {
		{"elastic", {"energy","angle"}},
		{"inelastic", {"energy","omega0"}}
	};
	std::map<std::string,std::function<double(double)>> m_type_differential_map = {
		{"elastic", [](double theta){return 2.0*constant::pi*std::sin(theta);}},
		{"inelastic", [](double x){return 1.0;}}
	};
};

#endif
