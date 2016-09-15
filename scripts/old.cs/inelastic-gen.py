#!/usr/bin/env python3
# Based on Ashley/KieftBosch
# See /doc/?
# Use on KieftBosch optical ELF data (data/elf/*.dat)
import argparse
import numpy as np
from numpy import (log, sqrt, log10, pi)
import scipy.interpolate
from constants import (q_e, a_0)


def getdata(elf_file):
    w0_ls = []
    elf_ls = []
    with open(elf_file) as f:
        for line in f:
            vals = line.split()
            if len(vals) != 2:
                continue

            w0 = float(vals[0])     # in eV
            elf = float(vals[1])    # in ?

            if w0 < 0 or elf < 0:
                continue

            w0_ls.append(w0)
            elf_ls.append(elf)

    return np.array(w0_ls), np.array(elf_ls)


def L_Kieft(w0, K, F):
    a = w0 / K
    x1 = 2/a * (1 + sqrt(1 - 2*a)) - 1
    x2 = (K - F - w0) / (K - F + w0)
    L1 = log(x1) + log(x2)
    L2 = -log(a)
    return (w0 / K < 1/2) * 1.5*L1 + (w0 > 50) * L2


def L_Ashley_w_ex(w0, K, _):
    a = w0/K
    return (1 - a)*log(4/a) - 7/4*a + a**(3/2) - 33/32*a**2


def L_Ashley_wo_ex(w0, K, _):
    a = w0/K
    s = sqrt(1-2*a)
    return log((1 - a/2 + s)/(1 - a/2 - s))

methods = {
    'Kieft': L_Kieft,
    'Ashley_w_ex': L_Ashley_w_ex,
    'Ashley_wo_ex': L_Ashley_wo_ex
}

parser = argparse.ArgumentParser(
    description='Calculate inelastic cross-sections from optical data.')
parser.add_argument(
    '--elf-file', type=str, required=True,
    help='optical data energy loss file')
parser.add_argument(
    '--number-density', type=float, required=True,
    help='number density [#/m^3]')
parser.add_argument(
    '--fermi', type=float, default=0,
    help='Fermi energy [eV]')
parser.add_argument(
    '--L', choices=list(methods.keys()), default='Kieft',
    help='L')

args = parser.parse_args()

F = args.fermi
rho_n = args.number_density
L = methods[args.L]
w0, elf = getdata(args.elf_file)
elf_interp = scipy.interpolate.PchipInterpolator(w0, elf)

K_bounds = (max(2*F,1),1e6)

print('<cstable type="inelastic">')
for K in np.logspace(log10(K_bounds[0]), log10(K_bounds[1]), num=1024):
    print('<cross-section energy="{}*eV">'.format(K))
    for w in np.logspace(log10(w0[0]), log10(K/2 - 0.1), num=1024):
        dcs = elf_interp(w)*L(w, K, F) / (2*pi*K*a_0*q_e * rho_n)
        print('\t<insert omega0="{}*eV" dcs="{}*m^2" />'.format(w, dcs))
    print('</cross-section>')
print('</cstable>')
