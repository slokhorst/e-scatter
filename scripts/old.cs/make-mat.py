#!/usr/bin/env python3
import os
import sys


el_phon_fn = os.path.join(mat_dir, 'elastic-phonon.xml')
el_mott_fn = os.path.join(mat_dir, 'elastic-mott.xml')
el_fn = os.path.join(mat_dir, 'elastic.xml')
inel_fn = os.path.join(mat_dir, 'inelastic.xml')
ion_fn = os.path.join(mat_dir, 'ionization.xml')
osi_fn = os.path.join(mat_dir, 'outer_shell.dat')

for elem in mat['elements']:
    # Z = mat['elements'][elem]['Z']
    # count = mat['elements'][elem]['count']

    # element_el_mott_fn = os.path.join(
        # mat_dir, 'elastic-mott-{}.xml'.format(elem))
    # print('running Elsepa for element {}'.format(elem), file=sys.stderr)
    # run_elsepa(Z, element_el_mott_fn)
    # with open('tmp.xml', 'w') as tmp_xml:
        # subprocess.run(
            # ['bin/cstool', 'mad', str(1), el_mott_fn, str(count),
             # element_el_mott_fn],
             # stdout=tmp_xml, check=True)
    # os.rename('tmp.xml', el_mott_fn)

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
