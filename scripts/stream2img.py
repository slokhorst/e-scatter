#!/usr/bin/env python3
import sys
import stream
import numpy
import matplotlib.image
import matplotlib.cm
from collections import defaultdict

def parse(filename):
    se_map = defaultdict(int)
    bse_map = defaultdict(int)
    for particle in stream.read_all(filename):
        if particle.K > 50:
            bse_map[(particle.xi, particle.yi)] += 1
        else:
            se_map[(particle.xi, particle.yi)] += 1
    return se_map, bse_map


def to_matrix(i_map):
    sx, sy = 1, 1
    max_I = 1
    for (xi, yi), I in i_map.items():
        sx = max(sx, xi)
        sy = max(sy, yi)
        max_I = max(max_I, I)
    img = numpy.ndarray((sy+1, sx+1))
    for (xi, yi), I in i_map.items():
        img[yi, xi] = I/max_I
    return img

if len(sys.argv) != 3:
    print("usage: {} <input> <output_base>".format(sys.argv[0]))
    sys.exit(1)

se_map, bse_map = parse(sys.argv[1])
basename = sys.argv[2]
se_mat = to_matrix(se_map)
bse_mat = to_matrix(bse_map)

matplotlib.image.imsave(
    "{}_se.png".format(basename), se_mat, cmap=matplotlib.cm.gray)
matplotlib.image.imsave(
    "{}_bse.png".format(basename), bse_mat, cmap=matplotlib.cm.gray)