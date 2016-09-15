from noodles.run.run_with_prov import run_parallel_opt
from noodles.display import NCDisplay
from noodles.serial import Registry, Serialiser, base
from noodles.serial.numpy import registry as numpy_registry

from cstool.parse_input import read_input, pprint_settings, cstool_model
from cstool.mott import s_mott_cs

from cslib.dataframe import DataFrame, DCS
from cslib import units

import numpy as np


class SerQuantity(Serialiser):
    def __init__(self):
        super(SerQuantity, self).__init__('<quantity>')

    def encode(self, obj, make_rec):
        return make_rec(obj.to_tuple())

    def decode(self, cls, data):
        return units.Quantity.from_tuple(data)


class SerUnit(Serialiser):
    def encode(self, obj, make_rec):
        return make_rec(str(obj))

    def decode(self, cls, data):
        return units(data)


class SerStandardObject(Serialiser):
    def __init__(self, cls, items):
        super(SerStandardObject, self).__init__(cls)
        self.items = items

    def encode(self, obj, make_rec):
        return make_rec({k: getattr(obj, k) for k in self.items})

    def decode(self, cls, data):
        return cls(**data)


def quantity_hook(obj):
    if isinstance(obj, units.Quantity):
        return '<quantity>'

    if type(obj).__name__ == 'Quantity':
        return '<quantity>'

    if isinstance(obj, units.Unit):
        return '<unit>'

    return None


def registry():
    return Registry(
        parent=base() + numpy_registry(),
        types={
            DataFrame: SerStandardObject(DataFrame, ['data', 'units', 'comments']),
            DCS: SerStandardObject(DCS, ['energy', 'angle', 'cs'])
        },
        hooks={
            '<quantity>': SerQuantity(),
            '<unit>': SerUnit('<unit>')
        },
        hook_fn=quantity_hook)


if __name__ == "__main__":
    s = read_input("./materials/silicon.json")

    print(pprint_settings(cstool_model, s))
    print()
    print("Phonon loss: {:~P}".format(s.phonon_loss))
    print("Total molar weight: {:~P}".format(s.M_tot))
    print("Number density: {:~P}".format(s.rho_n))
    print()
    print("Computing Mott cross-sections using ELSEPA.")

    e = np.logspace(1, 5, 145) * units.eV
    f_mcs = s_mott_cs(s, e)

    with NCDisplay() as display:
        mcs = run_parallel_opt(
            f_mcs, n_threads=4, registry=registry,
            jobdb_file='cache.json', display=display)

    mcs.save_gnuplot('{}_mott.bin'.format(s.name))
