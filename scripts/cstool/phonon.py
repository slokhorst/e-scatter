# Based on Schreiber & Fitting
# See /doc/extra/phonon-scattering.lyx

from cslib import units, Q_
from cslib.numeric import (interpolate, identity)

from cstool.parse_input import read_input

from math import pi

import numpy as np
import sys

from numpy import (cos, exp, log10)
from functools import partial


def phonon_crosssection(eps_ac, c_s, M, rho_m, lattice=None, E_BZ=None):
    if lattice is None and E_BZ is None:
        raise ValueError("One of `lattice` and `E_BZ` should be given.")

    T_room = 300 * units.K
    E_BZ = E_BZ or units.h**2 / (2*units.m_e * lattice**2)

    A = 5*E_BZ
    rho_n = Q_('N_A') / M * rho_m
    # h_bar_w_BZ = Q_('h') * c_s / lattice
    # n_BZ = 1 / (exp(h_bar_w_BZ / units.k / T_room) - 1)
    sigma_ac = (Q_('m_e² k') * eps_ac**2 * T_room) / \
        (units.hbar**4 * c_s**2 * rho_m * rho_n)

    def dcs_lo(theta, E):
        """Phonon cross-section for low energies."""
        return 1

    def dcs_hi(theta, E):
        """Phonon cross-section for high energies.

        :param E: energy in Joules.
        :param theta: angle in radians."""
        return (4*A / E_BZ) * (1 - cos(theta))/2 * E/A

    def norm(theta, E):
        return sigma_ac / (4*pi) / (1 + (1 - cos(theta))/2 * E/A)**2

    def h(x):
        return -4*x**3 + 6*x**2 - 1

    def dcs(theta, E):
        g = interpolate(
            partial(dcs_lo, theta), partial(dcs_hi, theta),
            h, E_BZ / 4, E_BZ)
        return g(E) * norm(theta, E)

    # should have units of m²/sr
    return dcs


if __name__ == "__main__":
    import argparse
    import io

    parser = argparse.ArgumentParser(
        description='Calculate elastic phonon cross-sections for a material.')
    parser.add_argument(
        'material_file', type=str,
        help="Filename of material in JSON format.")
    args = parser.parse_args()

    s = read_input(args.material_file)
    if 'M_tot' not in s:
        s.M_tot = sum(e.M * e.count for e in s.elements.values())

    E_range = np.logspace(log10(0.01), log10(1000), num=100) * units.eV
    theta_range = np.linspace(0, pi, num=100) * units.rad

    csf = phonon_crosssection(s.eps_ac, s.c_s, s.M_tot, s.rho_m, s.lattice)
    cs = csf(theta_range[:, None], E_range)

    gp_bin = np.zeros(dtype='float32', shape=[101, 101])
    gp_bin[0, 0] = 100
    gp_bin[1:, 0] = theta_range
    gp_bin[0, 1:] = E_range
    gp_bin[1:, 1:] = cs
    gp_bin.tofile('test.bin')
