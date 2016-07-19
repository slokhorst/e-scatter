#!/usr/bin/python3
import sys
import struct
import numpy
import matplotlib.image
import matplotlib.cm
from collections import defaultdict


def parse(filename):
    se_map = defaultdict(int)
    bse_map = defaultdict(int)
    with open(filename, "rb") as f:
        while f.read(1):
            f.seek(-1, 1)
            x, y, z, dx, dy, dz, K = struct.unpack('7f', f.read(7*4))
            xi, yi = struct.unpack('2i', f.read(2*4))
            if K > 50:
                bse_map[(xi, yi)] += 1
            else:
                se_map[(xi, yi)] += 1
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

se_mat_l, se_mat_r = numpy.split(se_mat, 2, axis=1)
bse_mat_l, bse_mat_r = numpy.split(bse_mat, 2, axis=1)

matplotlib.image.imsave(
    "{}_se_l.png".format(basename), se_mat_l, cmap=matplotlib.cm.gray)
matplotlib.image.imsave(
    "{}_se_r.png".format(basename), se_mat_r, cmap=matplotlib.cm.gray)
matplotlib.image.imsave(
    "{}_bse_l.png".format(basename), bse_mat_l, cmap=matplotlib.cm.gray)
matplotlib.image.imsave(
    "{}_bse_r.png".format(basename), bse_mat_r, cmap=matplotlib.cm.gray)
