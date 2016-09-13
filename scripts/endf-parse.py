#!/usr/bin/env python3
# Use on
# EADL (Evaluated Atomic Data Library) and
# EEDL (Evaluated Electron Data Library) files
# in ENDF (Evaluated Nuclear Data File) format
# http://www.nndc.bnl.gov/endf/
# ENDF/B-VII.1 (2012):
# http://www.nndc.bnl.gov/endf/b7.1/download.html
import os
import sys
import glob

from cslib import units


def parse_endf_line(line):
    """Takes a line from an ENDF file and splits it in parts. The
    badly represented numbers are reformatted to comply with standards.
    In this process it is assumed that there are no negative numbers
    in the ENDF file."""
    return {
        0: line[0:11].replace('+', 'E+').replace('-', 'E-'),
        1: line[11:22].replace('+', 'E+').replace('-', 'E-'),
        2: line[22:33].replace('+', 'E+').replace('-', 'E-'),
        3: line[33:44].replace('+', 'E+').replace('-', 'E-'),
        4: line[44:55].replace('+', 'E+').replace('-', 'E-'),
        5: line[55:66].replace('+', 'E+').replace('-', 'E-'),
        "mat": int(line[66:70]),
        "mf": int(line[70:72]),
        "mt": int(line[72:75]),
        "ln": int(line[75:81])
    }


def parse_eadl_file(filename):
    """Parse a single EADL file, extracts only the binding energy
    and occupancy of each shell."""
    shells = {}

    with open(filename) as f:
        lines = f.readlines()
        ln = 0
        while ln < len(lines):
            row = parse_endf_line(lines[ln])
            # 28/533: Atomic Relaxation Data for Electrons and Photons
            if row["mf"] == 28 and row["mt"] == 533:
                nss = int(row[4])   # number of subshells
                ln += 1

                for i_s in range(nss):
                    row1 = parse_endf_line(lines[ln+0])
                    row2 = parse_endf_line(lines[ln+1])

                    subi = int(float(row1[0]))       # subshell designator
                    ebi = float(row2[0]) * units.eV  # binding energy (eV)
                    eln = float(row2[1])             # occupancy
                    shells[subi] = (ebi, eln)

                    ln += int(float(row1[5]))+2
            ln += 1

    return shells


def parse_eedl_file(filename):
    """Parse EEDL file, extracts the ionization cross-sections."""
    shell_cross_sections = {}

    with open(filename) as f:
        lines = f.readlines()
        ln = 0
        while ln < len(lines):
            row = parse_endf_line(lines[ln])
            # 23/534--23/599: Electroionization Subshell Cross
            # Sections shell 1--66
            if row["mf"] == 23 and (row["mt"] >= 534 and row["mt"] <= 599):
                row2 = parse_endf_line(lines[ln+1])
                # row3 = parse_endf_line(lines[ln+2])
                subi = row["mt"]-533
                # ebi = float(row2[0])    # binding energy (eV)
                n = int(row2[5])
                ln += 2

                col = 3
                for _ in range(n):
                    if col > 2:
                        ln += 1
                        rowd = parse_endf_line(lines[ln])
                        col = 0

                    energy = float(rowd[2*col+0]) * units.eV  # energy
                    cs = float(rowd[2*col+1]) * units.barn    # crosssection

                    if subi not in shell_cross_sections:
                        shell_cross_sections[subi] = {}
                    shell_cross_sections[subi] = (energy, cs)

                    col += 1
            ln += 1

    return shell_cross_sections


def parse_folder(dirname, z):
    eadl_filenames = os.path.join(
        dirname, "atomic_relax", "atom-{0:03d}_*_*.endf".format(z))
    eadl_files = glob.glob(eadl_filenames)
    eedl_filenames = os.path.join(
        dirname, "electrons", "e-{0:03d}_*_*.endf".format(z))
    eedl_files = glob.glob(eedl_filenames)

    if len(eadl_files) != 1:
        raise FileNotFoundError(eadl_filenames)
    if len(eedl_files) != 1:
        raise FileNotFoundError(eedl_filenames)

    eadl_fpath = eadl_files[0]
    eedl_fpath = eedl_files[0]

    shell_cross_sections = parse_eedl_file(eedl_fpath)

    for subi, (ebi, eln) in sorted(parse_eadl_file(eadl_fpath).items()):
        print('<cstable type="ionization" shell="{}" occupancy="{}"'
              ' binding-energy="{}*eV">'.format(subi, eln, ebi))
        for energy, cs in sorted(shell_cross_sections[subi].items()):
            print('\t<cross-section energy="{}*eV" cs="{}*m^2" />'
                  .format((energy-ebi), cs))
        print('</cstable>\n')

###

if len(sys.argv) != 3:
    print("usage: {} <ENDF-root-dir> <Z-number>".format(sys.argv[0]))
    sys.exit(1)

parse_folder(sys.argv[1], int(sys.argv[2]))
