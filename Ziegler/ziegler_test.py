from npat import Ziegler, Isotope

if __name__=='__main__':

	######## Ziegler Tests #############

	## Compound must be specified, and enough info to determine areal density
	## Units:
	## thickness: mm
	## mass: g
	## area: cm^2
	## ad (areal density): mg/cm^2
	## density: g/cm^3
	import numpy as np
	names = ['-10', '-5', '0', '+5', '+10' ]
	scaling_parameter = [0.9, 0.95, 1.0, 1.05, 1.1]

	#SS_316 stainless steel
	for i in range(len(scaling_parameter)):
		stack = [{'compound':'SS_316', 'name':'SS01', 'ad':102.0},

				{'compound':'Ni', 'name':'Ni01', 'ad':22.772},
				{'compound':'Ir', 'name':'Ir01', 'ad':55.174},
				{'compound':'Cu', 'name':'Cu01', 'ad':22.338},

				{'compound':'Fe', 'name':'Fe01', 'ad':20.030*scaling_parameter[i]},


				{'compound':'Ni', 'name':'Ni02', 'ad':23.118},
				{'compound':'Ir', 'name':'Ir02', 'ad':55.601},
				{'compound':'Cu', 'name':'Cu02', 'ad':22.325},

				{'compound':'Fe', 'name':'Fe02', 'ad':20.017*scaling_parameter[i]},


				{'compound':'Ni', 'name':'Ni03', 'ad':22.338},
				{'compound':'Ir', 'name':'Ir03', 'ad':55.643},
				{'compound':'Cu', 'name':'Cu03', 'ad':22.313},

				{'compound':'Fe', 'name':'Fe03', 'ad':19.948*scaling_parameter[i]},

				{'compound':'Ni', 'name':'Ni04', 'ad':20.704},
				{'compound':'Ir', 'name':'Ir04', 'ad':56.000},
				{'compound':'Cu', 'name':'Cu04', 'ad':22.284},


				{'compound':'Ni', 'name':'Ni05', 'ad':21.768},
				{'compound':'Ir', 'name':'Ir05', 'ad':55.161},
				{'compound':'Cu', 'name':'Cu05', 'ad':22.443},


				{'compound':'Ni', 'name':'Ni06', 'ad':22.861},
				{'compound':'Ir', 'name':'Ir06', 'ad':55.731},
				{'compound':'Cu', 'name':'Cu06', 'ad':22.396},


				{'compound':'Ni', 'name':'Ni07', 'ad':23.092},
				{'compound':'Ir', 'name':'Ir07', 'ad':56.685},
				{'compound':'Cu', 'name':'Cu07', 'ad':22.320},


				{'compound':'Ni', 'name':'Ni08', 'ad':22.409},
				{'compound':'Ir', 'name':'Ir08', 'ad':58.030},
				{'compound':'Cu', 'name':'Cu08', 'ad':22.401},


				{'compound':'Ni', 'name':'Ni09', 'ad':21.741},
				{'compound':'Ir', 'name':'Ir09', 'ad':56.669},
				{'compound':'Cu', 'name':'Cu09', 'ad':22.425},


				{'compound':'Ni', 'name':'Ni10', 'ad':23.093},
				{'compound':'Ir', 'name':'Ir10', 'ad':55.065},
				{'compound':'Cu', 'name':'Cu10', 'ad':22.314},

				{'compound':'SS_316', 'name':'SS02', 'ad':102.0} ]

		## must specify stack and beam...
		## stack is list of dicts, which specify composition of each
		## element of the stack
		##
		## beam is dict, specifying the isotope and incident energy
		## of the incoming beam.  dE0 (1-sigma width) and N (number of particles)
		## can optionally be included
		zg = Ziegler(stack=stack, beam_istp='2H', E0=33.0, N=1E6, max_steps=100)

		## name from stack is used for plotting - either incl_no_names=True must
		## be set, or name must be specified
		## if zg.plot() is called, all foils in stack are plotted that have a name.
		## optionally the first argument to plot is a list of samples that will be
		## regex matched to the names in the stack

		#zg.plot(['Ir'])#,'Ni','Ti'])

		## summarize either prints out the mean and 1-sigma energies of each
		## foil in the stack (with a name), or saves this info to a .csv using
		## the saveas argument.
		#zg.summarize(['Ir'])

		#zg.saveas('E_foils_test3.csv')


		zg.saveas('E_foils_Fe_{}.csv'.format(names[i]))

	## The results of the stack calculation can be found in the zg.stack
	## variable, which is a list of Sample objects which store the energy,
	## flux, and other info about each sample in the stack calculation.
	#print(zg.stack[2].meta)
	#print(zg.stack[1].energy)
	#print(zg.stack[1].flux)


	####### Isotope Tests ############

	# ISTP = Isotope('133BA')
	# print ISTP.name
	# print ISTP.element
	# print ISTP.A
	# print ISTP.isomer
	# print ISTP.isotope
	# print ISTP.E_level
	# print ISTP.TeX
	# print ISTP.mass
	# print ISTP.half_life(ISTP.optimum_units(),unc=True),ISTP.optimum_units()
	# print ISTP.gammas()
	# print ISTP.electrons()
	# print ISTP.beta_minus()
	# print ISTP.beta_plus()
	# print ISTP.alphas()
