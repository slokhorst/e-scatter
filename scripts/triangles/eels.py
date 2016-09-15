#!/usr/bin/env python3
import argparse
import math
import stream
import subprocess
import numpy as np
import matplotlib.pyplot as plt

plot_range=(0,100) #eV
bin_size=1 #nm

def teels(material_file, K0, num_pri, theta_c):
	with open('primary.stream', 'w') as f:
		p = stream.record(0,0,1,0,0,-1,K0,1,1)
		for i in range(0,num_pri):
			p.write(f.buffer)

	with open('detected.stream', 'w') as f:
		subprocess.run(['bin/cdsem', '../data/teels.tri', 'primary.stream',
		                material_file,
		                '--generate-secondary=0'],
		                stdout=f)

	w0_energy_losses = []
	for particle in stream.read_all('detected.stream'):
		if(particle.dz > 0):
			raise Exception("upwards electron detected")
		r = math.sqrt(particle.dx*particle.dx + particle.dy*particle.dy)
		theta = math.atan(r/(-particle.dz))

		if (theta < theta_c):
			w0_energy_losses.append(K0-particle.K)

	hist, bins = np.histogram(w0_energy_losses,
		bins=round((plot_range[1]-plot_range[0])/bin_size),
		range=plot_range)
	width = (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	plt.bar(center, hist/num_pri, align='center', width=width)
	plt.title('TEELS: zero momentum energy loss function ($\\theta<{}$rad)'.format(theta_c))
	plt.xlabel('energy loss [eV]')
	plt.show()

def reels(material_file, K0, num_pri):
	with open('primary.stream', 'w') as f:
		p = stream.record(0,0,1,0,0,-1,K0,1,1)
		for i in range(0,num_pri):
			p.write(f.buffer)

	with open('detected.stream', 'w') as f:
		subprocess.run(['bin/cdsem', '../data/reels.tri', 'primary.stream',
		                material_file,
		                '--generate-secondary=1'],
		                stdout=f)

	energy_losses = []
	for particle in stream.read_all('detected.stream'):
		if(particle.dz < 0):
			raise Exception("downwards electron detected")
		energy_losses.append(K0-particle.K)

	hist, bins = np.histogram(energy_losses, bins=1000)
	width = (bins[1] - bins[0])
	center = (bins[:-1] + bins[1:]) / 2
	plt.bar(center, hist/num_pri, align='center', width=width)
	plt.title('REELS: energy loss function')
	plt.xlabel('energy loss [eV]')
	plt.show()


parser = argparse.ArgumentParser(
	description='Generate an electron energy loss spectrum')
parser.add_argument(
	'material', type=str, help='material file')
parser.add_argument(
    '--type', choices=list(['TEELS','REELS']), required=True,
    help='transmissive or reflective electron loss spectrum')
parser.add_argument(
    '-k', type=float, default=10000,
    help='primary electron energy [eV]')
parser.add_argument(
    '-n', type=float, default=1000000,
    help='number of primary electrons')
parser.add_argument(
    '--theta-c', type=float, default=0.1,
    help='collection angle for TEELS [rad]')
args = parser.parse_args()

if(args.type == 'TEELS'):
	teels(args.material, args.k, args.n, args.theta_c)
else:
	reels(args.material, args.k, args.n)
