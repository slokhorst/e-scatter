{
	'name': 'pmma',

	'rho_m': 1.192e3,
	'fermi': 0*eV,
	'work_func': 2.5*eV,
	'band_gap': 5.6*eV,
	'lattice': 0.543e-9,  # silicon
	'c_s': 9040,          # silicon
	'eps_ac': 9.2*eV,     # silicon

	'inelastic_xml_file': 'inelastic-pmma.xml',

	'elements': {
		'H':  {'count': 8, 'Z': 1, 'M':  1.008e-3},
		'C':  {'count': 5, 'Z': 6, 'M': 12.011e-3},
		'O':  {'count': 2, 'Z': 8, 'M': 15.999e-3}
	}
}
