from collections import OrderedDict
import io
import numpy as np

from . import units


class DataFrame(object):
    """A class like the Pandas DataFrame; this one supports physical units
    using the `pint` module, and storage to HDF5 files.

    This is a wrapper around a higher dimensional Numpy array.

    The class supports both column and
    row based access; however, it is optimised to handle entire columns
    of data more efficiently.

    Rows start counting at 0. Every column must have a unit."""
    def __init__(self, data, units=None, comments=None):
        self.data = data
        self.units = units
        self.comments = comments
        self.unit_dict = OrderedDict(zip(self.data.dtype.names, units))

    def __getitem__(self, x):
        s = self.data[x]
        if isinstance(x, str):
            return s * self.unit_dict[x]
        elif isinstance(x, tuple) and isinstance(x[1], int):
            return s * self.units[x[1]]
        elif isinstance(x, tuple) and isinstance(x[1], slice):
            return DataFrame(s, self.units[x[1]])

    def __len__(self):
        return len(self.data)

    def __str__(self):
        of = io.BytesIO()
        np.savetxt(of, self.data, fmt='%.4e')
        return '# ' + ', '.join(
            '{0} ({1:~})'.format(n, u)
            for n, u in self.unit_dict.items()) + \
            '\n' + of.getvalue().decode()


class DCS(object):
    """Differential cross-section: energy, angle, cs.

    The energy and angle are 1d array quantities, with dimensions
    `N` and `M` respectively and having units with dimensionality
    of energy and angle.

    The cross-section is a 2d array of shape [N, M], having dimensionality
    of area."""
    def __init__(self, energy, angle, cs):
        self.energy = energy
        self.angle = angle
        self.cs = cs

        self._E_log_steps = np.log(self.energy[1:]/self.energy[:-1])
        self._a_steps = self.angle[1:] - self.angle[:-1]

        assert energy.dimensionality == units.J.dimensionality, \
            "Energy units check."
        assert angle.dimensionality == units.rad.dimensionality, \
            "Angle units check."
        assert cs.dimensionality == (units.m**2).dimensionality, \
            "Cross-section units check."
        assert cs.shape == (energy.size, angle.size), \
            "Array dimensions do not match."

    def save_gnuplot(self, filename):
        xlen, ylen = self.cs.shape
        gp_bin = np.zeros(dtype='float32', shape=[xlen+1, ylen+1])
        gp_bin[0, 0] = xlen
        gp_bin[1:, 0] = self.energy.to('eV')
        gp_bin[0, 1:] = self.angle.to('rad')
        gp_bin[1:, 1:] = self.cs.to('cm²')
        gp_bin.transpose().tofile(filename)

    def __call__(self, E, a):
        """Multi-linear interpolation on this DCS table.
        The interpolation is log-linear on energy and linear on angle."""
        E_idx = np.searchsorted(self.energy, E) - 1
        E_w = np.log(E / self.energy[E_idx]) / self._E_log_steps[E_idx]

        a_idx = np.searchsorted(self.angle, a) - 1
        a_w = (a - self.angle[a_idx]) / self._a_steps[a_idx]

        idx = np.meshgrid(E_idx, a_idx).transpose([1, 2, 0])
        new_cs = (1 - E_w) * (1 - a_w) * self.cs[idx] \
            + E_w * (1 - a_w) * self.cs[idx + [1, 0]] \
            + (1 - E_w) * a_w * self.cs[idx + [0, 1]] \
            + E_w * a_w * self.cs[idx + [1, 1]]
        return DCS(E, a, new_cs)
