#!/usr/bin/python3
# Based on Schreiber & Fitting
# See /doc/extra/phonon-scattering.lyx
# WARNING: line 44, 64: added pi's to match Kieft
import argparse
import math
import numpy
import sys

h = 6.626070040e-34
h_bar = h/(2*math.pi)
k_B = 1.3806488e-23
m_e = 9.10938215e-31
q_e = 1.602176565e-19
eV = q_e
N_A = 6.022141e23
T=300

parser = argparse.ArgumentParser(description='Calculate elastic phonon cross-sections for a material.')
parser.add_argument('--eps_ac', type=float, required=True, help='acoustic deformation potential [eV]')
parser.add_argument('--c_s', type=float, required=True, help='speed of sound in material [m/s]')
parser.add_argument('--M', type=float, required=True, help='molar weight [kg/mol]')
parser.add_argument('--rho_m', type=float, required=True, help='mass density [kg/m^3]')
parser.add_argument('--a', type=float, help='lattice constant [m]')
parser.add_argument('--e_bz', type=float, help='brioullin zone energy [J] (can be deduced from a)')

args = parser.parse_args()
eps_ac = args.eps_ac
c_s = args.c_s
rho_m = args.rho_m
rho_n = N_A/args.M*rho_m
if args.e_bz is not None:
	if args.a is not None:
		print("WARNING: ignoring parameter a", file=sys.stderr)
	E_BZ = args.e_bz
elif args.a is not None:
	E_BZ = h**2/(2*m_e*args.a**2)
else:
	raise SyntaxError("must define either a or e_bz")

A = 5*E_BZ

sigma_ac = (m_e**2*eps_ac**2*k_B*T)/(math.pi*h_bar**4*c_s**2*rho_m*rho_n)
sigma_ac = sigma_ac*math.pi # Kieft?! normalization

def dcs(E,theta):
	if E < E_BZ/4:
		return dcs_lo(E,theta)
	elif E > E_BZ:
		return dcs_hi(E,theta)
	else:
		x = (E-E_BZ/4)/(E_BZ-E_BZ/4)
		return (1-x)*dcs_lo(E_BZ/4,theta) + x*dcs_hi(E_BZ,theta)


def dcs_lo(E,theta):
	return sigma_ac/(4*math.pi) \
		* 1/(1+(1-math.cos(theta))/2*E/A)**2


def dcs_hi(E,theta):
	return sigma_ac/(4*math.pi) \
		* (4*A)/(E_BZ) * ((1-math.cos(theta))/2*E/A)/(1+(1-math.cos(theta))/2*E/A)**2 \
		* math.pi # Kieft!? magically appeared in their derivation

###

print('<cstable type="elastic">')
for E in numpy.logspace(math.log10(0.01*eV), math.log10(1000*eV), num=100):
	print('\t<cross-section energy="{energy}*eV">'.format(energy=E/eV))
	for theta in numpy.linspace(0, math.pi, num=100):
		print('\t\t<insert angle="{angle}" dcs="{dcs}*m^2/sr" />'.format(angle=theta,dcs=dcs(E,theta)))
	print('\t</cross-section>')
print('</cstable>\n')
