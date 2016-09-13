from cslib import units
from cstool.parse_input import cstool_model, read_input, pprint_settings


class ObjectUnion(object):
    def __init__(self, *args):
        self.s = args

    def __getitem__(self, key):
        for d in self.s:
            if key in dir(d):
                value = getattr(d, key)
                return value
        else:
            raise AttributeError("{} not in union of dicts.".format(key))


if __name__ == "__main__":
    s = read_input("silicondioxide.json")

    print(eval("sum(e.M * e.count for e in elements.values())", {}, s))
    if 'M_tot' not in s:
        s.M_tot = sum(e.M * e.count for e in s.elements.values())

    print(eval("N_A / M_tot * rho_m", {}, ObjectUnion(s, units)))
    if 'rho_n' not in s:
        s.rho_n = (units.N_A / s.M_tot * s.rho_m).to("cm⁻³")

    print(eval("h * c_s / lattice", {}, ObjectUnion(s, units)))
    if 'phonon_loss' not in s:
        s.phonon_loss = (units.h * s.c_s / s.lattice).to("eV")

    s = pprint_settings(cstool_model, s)
    print(s)

