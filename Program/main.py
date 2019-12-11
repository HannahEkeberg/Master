
import numpy as np

A0 = 9003.74
mass_density = 22.772e20
lamb = 5.7664e-5
t = 3600
integral = 1.32e-26
ec = 1.60e-10

I = (A0*ec)/ (mass_density*(1-np.exp(-lamb*t))*integral)
print(I)



"""
#from single_decay_A0 import *
from foil_info import *
from beam_current_FoilReact import *
#from beam_current_nov19 import *
from beam_current_test_des19 import *
from ziegler_sorting import *

ziegler_file_SS_n10 = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_-10_fluxes.csv'
ziegler_file_SS_n5  = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_-5_fluxes.csv'
ziegler_file_SS_0   = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_0_fluxes.csv'
ziegler_file_SS_p5   = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_+5_fluxes.csv'
ziegler_file_SS_p10  = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_SS_+10_fluxes.csv'
filenames = [ziegler_file_SS_n10, ziegler_file_SS_n5,ziegler_file_SS_0, ziegler_file_SS_p5,ziegler_file_SS_p10]
names = ['-SS-10%','-SS-5%', '-SS0%', '-SS+5%', '-SS+10%']


def Beam_current_minimization(filenames, names, foil):

    for i in filenames:
        E_Ni, F_Ni, E_Cu, F_Cu, E_Fe, F_Fe = sort_ziegler(i)
        #name = 'Beam-Current-SS+10%'
        plot(names[i])
"""

#I =  linear_fit()[-1].T

#I = np.argsort(I[0])
#print(I)
