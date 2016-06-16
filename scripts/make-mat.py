#!/usr/bin/python3
# WARNING: line 103: Elastic cross-section from Elsepa multiplied by 2 to match Czyzewski
import os
import glob
import subprocess
import numpy
import math
import sys

N_A = 6.022141e23
q_e = 1.602176565e-19
eV = q_e

if len(sys.argv) != 2:
	print("usage: {} <material_definition.py>".format(sys.argv[0]))
	sys.exit(1)

mat_def_fn = sys.argv[1]
with open(mat_def_fn) as f:
	mat_def = f.read()
mat = eval(mat_def)

mat['M_tot'] = 0
for elem in mat['elements']:
	mat['M_tot'] += mat['elements'][elem]['count']*mat['elements'][elem]['M']
mat['rho_n'] = N_A/mat['M_tot'] * mat['rho_m']

###

work_dir = os.getcwd()
elsepa_dir = os.path.join(work_dir,"data/elsepa")
endf_dir = os.path.join(work_dir,"data/ENDF-B-VII.1")
mat_dir = os.path.join(work_dir,"data/materials",mat['name'])

###

def run_elsepa(Z, out_fn):
	out_xml = open(out_fn, 'w')
	os.chdir(elsepa_dir)
	for f in glob.glob('*.dat'):
		os.remove(f)

	no_muffin_Z = [1, 7, 8]

	muffin = 1
	if Z in no_muffin_Z:
		muffin = 0

	elscata_in = ''
	elscata_in += 'IZ      {}\n'.format(Z)      # atomic number                             [none]
	elscata_in += 'MNUCL   {}\n'.format(3)      # rho_n (1=P, 2=U, 3=F, 4=Uu)                [  3]
	elscata_in += 'MELEC   {}\n'.format(4)      # rho_e (1=TFM, 2=TFD, 3=DHFS, 4=DF, 5=file) [  4]
	elscata_in += 'MUFFIN  {}\n'.format(muffin) # 0=free atom, 1=muffin-tin model            [  0]
	elscata_in += 'IELEC   {}\n'.format(-1)     # -1=electron, +1=positron                   [ -1]
	elscata_in += 'MEXCH   {}\n'.format(1)      # V_ex (0=none, 1=FM, 2=TF, 3=RT)            [  1]
	elscata_in += 'MCPOL   {}\n'.format(0)      # V_cp (0=none, 1=B, 2=LDA)                  [  0]
	elscata_in += 'IHEF    {}\n'.format(0)      # high-E factorization (0=no, 1=yes, 2=Born) [  1]
	for E in [10, 20, 30, 40, 50, 60, 70, 80, 90, \
	          100, 200, 300, 400, 500, 600, 700, 800, 900, \
	          1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, \
	          10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000]:
		elscata_in += 'EV      {}\n'.format(E)  # kinetic energy (eV)                       [none]

	subprocess.run([os.path.join(elsepa_dir,'elscata')], input=bytes(elscata_in, 'UTF-8'), stdout=subprocess.DEVNULL)

	os.chdir(work_dir)

	subprocess.run(["scripts/elsepa-parse.py", elsepa_dir], stdout=out_xml)

def run_phonon(eps_ac, c_s, rho_m, M_tot, lattice, out_fn):
	out_xml = open(out_fn, 'w')
	subprocess.run(['scripts/phonon-gen.py','--eps_ac',str(eps_ac),'--c_s',str(c_s),'--rho_m',str(rho_m),'--M',str(M_tot),'--a',str(lattice)], stdout=out_xml)

def run_endf(Z, out_fn):
	out_xml = open(out_fn, 'w')
	subprocess.run(["scripts/endf-parse.py", endf_dir, str(Z)], stdout=out_xml)

###

el_phon_fn = os.path.join(mat_dir,'elastic-phonon.xml')
el_mott_fn = os.path.join(mat_dir,'elastic-mott.xml')
el_fn      = os.path.join(mat_dir,'elastic.xml')
inel_fn    = os.path.join(mat_dir,'inelastic.xml')
ion_fn     = os.path.join(mat_dir,'ionization.xml')

el_xml     = open(el_fn,'w')
inel_xml   = open(inel_fn,'w')
ion_xml    = open(ion_fn,'w')

with open(el_mott_fn,'w') as el_mott_xml:
	el_mott_xml.write('<cstable type="elastic">\n</cstable>\n');

i=0
for elem in mat['elements']:
	i+=1
	Z = mat['elements'][elem]['Z']
	count = mat['elements'][elem]['count']

	element_el_mott_fn = os.path.join(mat_dir,"elastic-mott-{}.xml".format(elem))
	run_elsepa(Z, element_el_mott_fn)
	with open('tmp.xml','w') as tmp_xml:
		subprocess.run(['bin/cstool','mad',str(1),el_mott_fn,str(count*2),element_el_mott_fn], stdout=tmp_xml) # WARNING 2 ADDED TO MATCH KIEFT
	os.rename('tmp.xml',el_mott_fn)

	element_ion_fn = os.path.join(mat_dir,"ionization-{}.xml".format(elem))
	run_endf(Z, element_ion_fn)
	subprocess.run(['bin/cstool','mad',str(count),element_ion_fn], stdout=ion_xml)
#subprocess.run(['scripts/inelastic-gen.py','--number-density',str(mat['rho_n']),'--elf-file',os.path.join(mat_dir,"elf.dat")], stdout=inel_xml)
subprocess.run(['cp',mat['inelastic_xml_file'],inel_fn])
run_phonon(mat['eps_ac'], mat['c_s'], mat['rho_m'], mat['M_tot'], mat['lattice'], el_phon_fn)
subprocess.run(['bin/cstool','merge',el_phon_fn,el_mott_fn,'100*eV','200*eV'], stdout=el_xml)

subprocess.run(['bin/compile-mat', \
	'{}/{}.mat'.format(mat_dir, mat['name']), \
	'--name', mat['name'], \
	'--elastic', el_fn, \
	'--inelastic', inel_fn, \
	'--ionization', ion_fn, \
	'--number-density', str(mat['rho_n']), \
	'--fermi-energy', str(mat['fermi']), \
	'--work-function', str(mat['work_func']), \
	'--band-gap', str(mat['band_gap']) \
])
