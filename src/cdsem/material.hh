/**
 * @file src/cdsem/material.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__MATERIAL__HEADER_INCLUDED
#define eSCATTER__CDSEM__MATERIAL__HEADER_INCLUDED

#include <map>
#include <string>
#include <vector>
#include <common/archive.hh>
#include <common/optional.hh>

class material {
friend archive::ostream& operator<<(archive::ostream&, const material&);
friend archive::istream& operator>>(archive::istream&, material&);
public:
    material() = default;
    material(const std::string& name, double fermi, double barrier, double density);
    /*!
     * @param[in] name
     *  Human readable identification string for the material.
     * @param[in] fermi
     *  Fermi energy.
     * @param[in] barrier
     *  Minimum energy required for an electron to escape from the material.
     * @param[in] bandgap
     *  Energy gap between the valence band and the conduction band.
     * @param[in] density
     *  The number density of the material.
     */
    material(const std::string& name, double fermi, double barrier, double bandgap, double density);
    inline const std::string& name() const;
    inline double fermi() const;
    inline double barrier() const;
    inline const optional<double>& bandgap() const;
    inline double density() const;
    inline double elastic_tcs(double K) const;
    inline double elastic_dcs(double K, double P) const;
    inline double inelastic_tcs(double K) const;
    inline double inelastic_dcs(double K, double P) const;
    double ionization_energy(double K, double P) const;
    material& set_elastic_data(double K, const std::map<double,double>& dcs_map);
    material& set_inelastic_data(double K, const std::map<double,double>& dcs_map);
    material& set_ionization_data(double B, const std::map<double,double>& tcs_map);
public:
    std::string _name;
    double _fermi = 0;
    double _barrier = 0;
    optional<double> _bandgap;
    double _density = 0;
    std::map<double,double> _elastic_tcs;
    std::map<double,std::map<double,double>> _elastic_dcs;
    std::map<double,double> _inelastic_tcs;
    std::map<double,std::map<double,double>> _inelastic_dcs;
    std::map<double,std::map<double,double>> _ionization_tcs;
};

archive::ostream& operator<<(archive::ostream&, const material&);
archive::istream& operator>>(archive::istream&, material&);

#include "material.inc"

#endif