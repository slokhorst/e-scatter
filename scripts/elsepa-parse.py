#!/usr/bin/python3
# Use on ELSEPA output
import os
import sys

def parse_tcstable_dat(filename):
	energy_tcs = {}

	with open(filename) as f:
		for line in f:
			if line.strip()[0] == "#":
				continue
			row = line.split()

			energy = float(row[0]) # energy (eV)
			tcs_cm2 = float(row[1]) # cross section (cm^2)
			tcs = tcs_cm2*1e-4

			energy_tcs[energy] = tcs

	return energy_tcs

def parse_dcs_dat(filename):
	angle_dcs = {}

	with open(filename) as f:
		for line in f:
			if line.strip()[0] == "#":
				continue
			row = line.split()

			angle = float(row[0]) # deflection angle (degrees)
			dcs_cm2 = float(row[2]) # differential cross section (cm^2/sr)
			dcs = dcs_cm2*1e-4

			angle_dcs[angle] = dcs

	return angle_dcs

def parse_folder(dirname):
	if not os.path.isdir(dirname):
		raise NotADirectoryError(dirname)

	tcstable_fpath = os.path.join(dirname,"tcstable.dat")
	if not os.path.isfile(tcstable_fpath):
		raise FileNotFoundError(tcstable_fpath)

	energy_tcs = parse_tcstable_dat(tcstable_fpath)

	print('<cstable type="elastic">')
	for energy, tcs in sorted(energy_tcs.items()):
		e_str = "{:8.3e}".format(energy).replace(".","p").replace("+","")
		dcs_fpath = os.path.join(dirname,"dcs_{}.dat".format(e_str))
		if not os.path.isfile(dcs_fpath):
			raise FileNotFoundError(dcs_fpath)

		print('\t<cross-section energy="{}*eV"><!--tcs="{}*m^2"-->'.format(energy,tcs))
		for angle, dcs in sorted(parse_dcs_dat(dcs_fpath).items()):
			print('\t\t<insert angle="{}*pi/180" dcs="{}*m^2/sr" />'.format(angle,dcs))
		print('\t</cross-section>')
	print('</cstable>\n')

###

if len(sys.argv) != 2:
	print("usage: {} <input-directory>".format(sys.argv[0]))
	sys.exit(1)

parse_folder(sys.argv[1])
