#!/usr/bin/python3
# From Ashley
# See /doc/?
import argparse
import math
import numpy
import sys

q_e = 1.602176565e-19
eV = q_e

parser = argparse.ArgumentParser(description='Calculate inelastic cross-sections from optical data.')
parser.add_argument('--elf-file', type=str, required=True, help='optical data energy loss file')
parser.add_argument('--number-density', type=float, required=True, help='numer density [#/m^3]')

args = parser.parse_args()

def dcs(E,omega0):
	return 0

###

print('<cstable type="inlelastic">')
for E in numpy.logspace(math.log10(0.01*eV), math.log10(1000*eV), num=100):
	omega0_max = 1
	print('\t<cross-section energy="{energy}*eV">'.format(energy=E/eV))
	for omega0 in numpy.linspace(0, omega0_max, num=100):
		print('\t\t<insert omega0="{omega0}" dcs="{dcs}*m^2" />'.format(omega0=omega0,dcs=dcs(E,omega0)))
	print('\t</cross-section>')
print('</cstable>\n')
