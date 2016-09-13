#!/usr/bin/env python3
import os
import glob
import subprocess
import sys
from constants import (N_A, eV, h)

# this script is being unified in the Script Unification Programme,
# all functionality that has been unified is signed
#   | SUP DONE  -> functionality is available in new framework
#   | SUP RED   -> functionality has become redundant

# command line | SUP RED
if len(sys.argv) != 2:
    print("usage: {} <material_definition.py>".format(sys.argv[0]))
    sys.exit(1)

# read settins from python file | SUP DONE
mat_def_fn = sys.argv[1]
with open(mat_def_fn) as f:
    mat_def = f.read()
mat = eval(mat_def)

# supply additional settings | SUP DONE
mat['M_tot'] = 0
for elem in mat['elements']:
    mat['M_tot'] += mat['elements'][elem]['count']*mat['elements'][elem]['M']
mat['rho_n'] = N_A/mat['M_tot'] * mat['rho_m']
mat['phonon_loss'] = h*mat['c_s']/mat['lattice']


work_dir = os.getcwd()
elsepa_dir = os.path.join(work_dir, "data/elsepa")
endf_dir = os.path.join(work_dir, "data/ENDF-B-VII.1")
mat_dir = os.path.join(work_dir, "data/materials", mat['name'])


# run ELSEPA | SUP DONE
def run_elsepa(Z, out_fn):
    pass


# phonon-gen has inconsistencies with the documentation
def run_phonon(eps_ac, c_s, rho_m, M_tot, lattice, out_fn):
    out_xml = open(out_fn, 'w')
    subprocess.run(
        ['scripts/phonon-gen.py', '--eps_ac', str(eps_ac), '--c_s', str(c_s),
         '--rho_m', str(rho_m), '--M', str(M_tot), '--a', str(lattice)],
        stdout=out_xml, check=True)


def run_endf(Z, out_fn):
    out_xml = open(out_fn, 'w')
    subprocess.run(
        ['scripts/endf-parse.py', endf_dir, str(Z)],
        stdout=out_xml, check=True)

###

el_phon_fn = os.path.join(mat_dir, 'elastic-phonon.xml')
el_mott_fn = os.path.join(mat_dir, 'elastic-mott.xml')
el_fn = os.path.join(mat_dir, 'elastic.xml')
inel_fn = os.path.join(mat_dir, 'inelastic.xml')
ion_fn = os.path.join(mat_dir, 'ionization.xml')
osi_fn = os.path.join(mat_dir, 'outer_shell.dat')

with open(el_mott_fn, 'w') as el_mott_xml:
    el_mott_xml.write('<cstable type="elastic">\n</cstable>\n')

for elem in mat['elements']:
    Z = mat['elements'][elem]['Z']
    count = mat['elements'][elem]['count']

    element_el_mott_fn = os.path.join(
        mat_dir, 'elastic-mott-{}.xml'.format(elem))
    print('running Elsepa for element {}'.format(elem), file=sys.stderr)
    run_elsepa(Z, element_el_mott_fn)
    with open('tmp.xml', 'w') as tmp_xml:
        subprocess.run(
            ['bin/cstool', 'mad', str(1), el_mott_fn, str(count),
             element_el_mott_fn],
             stdout=tmp_xml, check=True)
    os.rename('tmp.xml', el_mott_fn)

    element_ion_fn = os.path.join(mat_dir, 'ionization-{}.xml'.format(elem))
    print('parsing ENDF files for element {}'.format(elem), file=sys.stderr)
    run_endf(Z, element_ion_fn)
    with open(ion_fn, 'w') as ion_xml:
        subprocess.run(
            ['bin/cstool', 'mad', str(count), element_ion_fn],
            stdout=ion_xml, check=True)

print('generating phonon cross-sections', file=sys.stderr)
run_phonon(mat['eps_ac'], mat['c_s'], mat['rho_m'],
    mat['M_tot'], mat['lattice'], el_phon_fn)

with open(el_fn, 'w') as el_xml:
    subprocess.run(
        ['bin/cstool', 'merge', el_phon_fn, el_mott_fn, '100*eV', '200*eV'],
        stdout=el_xml, check=True)

with open('tmp.xml', 'w') as tmp_xml:
    subprocess.run(
        ['bin/cstool', 'shift', el_fn, str(mat['fermi'])],
        stdout=tmp_xml, check=True)
os.rename('tmp.xml', el_fn)

with open(mat['elf-file'], 'r') as elf_f, open(osi_fn, 'w') as osi_f:
    osi_f.write(elf_f.readline())

print('generating inelastic cross-sections', file=sys.stderr)
with open(inel_fn, 'w') as inel_xml:
    subprocess.run(
        ['scripts/inelastic-gen.py', '--number-density', str(mat['rho_n']),
         '--elf-file', mat['elf-file'], '--fermi', str(mat['fermi'] / eV)],
        stdout=inel_xml, check=True)

print('compiling material', file=sys.stderr)
subprocess.run(
    ['bin/cstool', 'compile-mat',
     '{}/{}.mat'.format(mat_dir, mat['name']),
     '--name', mat['name'],
     '--elastic', el_fn,
     '--inelastic', inel_fn,
     '--ionization', ion_fn,
     '--outer-shell', osi_fn,
     '--number-density', str(mat['rho_n']),
     '--fermi-energy', str(mat['fermi']),
     '--work-function', str(mat['work_func']),
     '--band-gap', str(mat['band_gap']),
     '--phonon-loss', str(mat['phonon_loss'])],
    check=True)
