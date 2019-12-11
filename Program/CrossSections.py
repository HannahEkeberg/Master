import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import N_A, elementary_charge
import pandas as pd
import os

from foil_info import *
from beam_current import *

I_weighted =125 #nA
irr_time = 3600 #seconds
lamb = Ni_56Ni()[-1]

def Fe():
    mol_mass_Fe = 55.8450 #g/mol
    mass_density = np.array((20.030, 20.017, 19.948))/1e3
    sigma_mass_density = np.array((0.110, 0.034, 0.114))/1e3
    mass_density = (mass_density*N_A)/mol_mass_Fe #nuclei/cm^2
    sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Fe
    return mass_density

def Cu():
    mol_mass_Cu = 63.546 #g/mol
    mass_density = np.array((22.338, 22.325, 22.313, 22.284, 22.443, 22.396, 22.320, 22.401, 22.425, 22.314)) /1e3
    sigma_mass_density = np.array((0.048, 0.028, 0.043, 0.027, 0.028, 0.012, 0.014, 0.033, 0.041, 0.047))/1e3
    mass_density = (mass_density*N_A)/mol_mass_Cu #nuclei/cm^2
    sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Cu
    return mass_density

def Ni():
    mol_mass_Ni = 58.69340 #g/mol
    mass_density = np.array((22.772, 23.118,22.338,20.704,21.768,22.861,23.092,22.409,21.741,23.093))/1e3   #g/cm^2, divide by 1e3 to get g instead of mg.
    sigma_mass_density = np.array((0.138, 0.096, 0.066, 0.068, 0.045, 0.123,0.078, 0.124, 0.073, 0.024))/1e3 #g/cm^2
    mass_density = (mass_density*N_A)/mol_mass_Ni #nuclei/cm^2
    sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Ni
    return mass_density


def cross_section(filename):
    mass_density = Ni()
    n = 10
    A0_array = np.zeros(n)
    sigma = np.zeros(n)
    for i in range(n):
        A0 = np.genfromtxt(filename, delimiter=',', usecols=[i])
        A0_array[i] = A0[0]
    for j in range(n):
        sigma[j] = A0_array[j] / (mass_density[j] * I_weighted*(1-np.exp(-lamb*irr_time)))   #I_weighted is suppose to be weighed chi squared for beam current.
    return sigma

def plot_CrossSections(x, y, reaction):
    plt.plot(x,y, '.')
    plt.plot(x,y, '--', linewidth=0.1)
    plt.xlabel('Beam energy (MeV)')
    plt.ylabel('Cross section (barns)')
    plt.title('Cross section for reaction {}'.format(reaction))
    plt.savefig('{}_CS'.format(reaction), dpi=300)
    plt.show()

#E_mean =
#E = weighted_average_beam_energy(E_mean)
E = np.linspace(0,33,10)
#CS = cross_section('Cu_65Zn.csv')
#CS = cross_section('Ir_193mPt.csv')
#CS = cross_section('Ni_56Ni.csv')
#CS = cross_section('Ni_58Co.csv')
#CS = cross_section('Ni_58mCo.csv')
#CS = cross_section('Ni_56Co.csv')
#CS = cross_section('Ni_62Zn.csv')
#CS = cross_section('Ni_61Cu.csv')
CS = cross_section('Ni_57Co.csv')
plot_CrossSections(E,CS, 'Ni_57Co')
