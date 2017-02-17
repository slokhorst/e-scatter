#!/usr/bin/env python3
import argparse
import sys
import numpy as np
from numpy import (random)
import struct

parser = argparse.ArgumentParser(
    description='Create exposure file containing primary electrons')
parser.add_argument(
    '--dose', '-D', type=float, default=2000.0,
    help='dose [uC/cm^2]')
parser.add_argument(
    '--spot', '-S', type=float, default=1.5,
    help='spot FW50 [nm]')
parser.add_argument(
    '--energy', '-K', type=float, default=300.0,
    help='energy [eV]')
parser.add_argument(
    '--nx', type=int, default=128,
    help='x resolution [px]')
parser.add_argument(
    '--ny', type=int, default=1024,
    help='y resolution [px]')
parser.add_argument(
    '--sx', type=float, default=64.0,
    help='x scan area (-sx/2 to +sx/2) [nm]')
parser.add_argument(
    '--sy', type=float, default=1000.0,
    help='y scan area (-sy/2 to +sy/2) [nm]')
parser.add_argument(
    '--pz', type=float, default=34.0,
    help='z position [nm]')
args = parser.parse_args()

X = np.linspace(-args.sx/2, +args.sx/2, num=args.nx)
Y = np.linspace(-args.sy/2, +args.sy/2, num=args.ny)
Z = args.pz
dx = 0.0
dy = 0.0
dz = -1.0

K = args.energy
D = args.dose*1e-6/(1e-2*1e-2)
FW50 = args.spot
std = FW50/(2*np.sqrt(2*np.log(2)))
n = (D*args.sx*1e-9*args.sy*1e-9)/1.60217662e-19
nperpx = n/(X.size*Y.size)

print('area: {} x {} @ z={} [nm]'.format(args.sx, args.sy, args.pz),
      file=sys.stderr)
print('resolution: {} x {} [px]'.format(args.nx, args.ny),
      file=sys.stderr)
print('spot size: {} [nm]'.format(args.spot),
      file=sys.stderr)
print('energy: {} [eV]'.format(args.energy),
      file=sys.stderr)
print('dose: {} [μC/cm^2]'.format(args.dose),
      file=sys.stderr)
print('      {} [e/px]'.format(nperpx),
      file=sys.stderr)

for yi in range(0, Y.size):
    print('generating exposure [{}%]'.format(int(100*yi/Y.size)),
          end="\r", file=sys.stderr)
    for xi in range(0, X.size):
        R = random.multivariate_normal(
            [X[xi], Y[yi]], std*std * np.identity(2),
            random.poisson(nperpx))
        for r in R:
            sys.stdout.buffer.write(
                struct.pack('7f', r[0], r[1], Z, dx, dy, dz, K))
            sys.stdout.buffer.write(
                struct.pack('2i', xi, yi))

sys.stdout.flush()
print("\n\r", end="", file=sys.stderr)
