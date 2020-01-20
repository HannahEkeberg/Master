import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import N_A, elementary_charge
import pandas as pd
import os

from foil_info import *
from beam_current_FoilReact import *
from ZieglerFiles import ziegler_files
from des19_BeamCurrent import *

files,names = ziegler_files()
##Ni-3% : 25
#print(names)
#myclass = BeamCurrent(files[25], sort_ziegler, Fe_foil, Ni_foil, Cu_foil)
#I = myclass.current_for_CS()

#print(Cu_foil('Cu_62Zn')[-1])
#path = os.getcwd() + '/activity_csv/'
#file = path + 'Ni_61Cu.csv'
#file = path + 'Ir_193mPt.csv'

#print(getA0(file, 10))

class CrossSections:

    def __init__(self,ziegler_file):
        self.ziegler_file = ziegler_file
        current_class = BeamCurrent(self.ziegler_file)
        #current_class = BeamCurrent(self.ziegler_file,)
        self.I = current_class.current_for_CS()
        t_irr = 3600
        self.I = np.true_divide(self.I, 3600.)   #?????

        #self.I = np.ones(10)*128.5
        self.E_Fe, self.E_Ni, self.E_Cu, self.E_Ir, self.dE_Fe, self.dE_Ni, self.dE_Cu, self.dE_Ir = current_class.current_for_CS(return_energies=True)
        self.irr_time       = 3600    #seconds
        self.sigma_irr_time = 3       #seconds
        self.path = os.getcwd() + '/activity_csv/'


    def get_variables(self, react_func, target_func):
        #react_func is eg. Fe_56Co()
        #target_func is eg. Ir_foil()
        lambda_ = react_func[-1]
        mass_density = target_func[-2]
        sigma_mass_density = target_func[-1]
        return lambda_, mass_density, sigma_mass_density

    def mass_density(self,foil):
        if foil == 'Ir':
            n=10
            mol_mass_Ir = 192.217 #g/mol
            mass_density = np.array((55.174, 55.601, 55.643, 56.000, 55.161, 55.731, 56.685, 58.030, 56.669, 55.065))/1e3
            sigma_mass_density = np.array((0.053, 0.238, 0.121, 0.109, 0.081, 0.088, 0.085, 0.130, 0.043, 0.055))/1e3
            mass_density = (mass_density*N_A)/mol_mass_Ir #nuclei/cm^2
            sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Ir
            E = self.E_Ir; dE = self.dE_Ir
        elif foil == 'Fe':
            n = 3 #number of foils
            mol_mass_Fe = 55.8450 #g/mol
            mass_density = np.array((20.030, 20.017, 19.948))/1e3
            sigma_mass_density = np.array((0.110, 0.034, 0.114))/1e3
            mass_density = (mass_density*N_A)/mol_mass_Fe #nuclei/cm^2
            sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Fe
            E = self.E_Fe; dE = self.dE_Fe
        elif foil == 'Ni':
            n = 10
            mol_mass_Ni = 58.69340 #g/mol
            mass_density = np.array((22.772, 23.118,22.338,20.704,21.768,22.861,23.092,22.409,21.741,23.093))/1e3   #g/cm^2, divide by 1e3 to get g instead of mg.
            sigma_mass_density = np.array((0.138, 0.096, 0.066, 0.068, 0.045, 0.123,0.078, 0.124, 0.073, 0.024))/1e3 #g/cm^2
            mass_density = (mass_density*N_A)/mol_mass_Ni #nuclei/cm^2
            sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Ni
            E = self.E_Ni; dE = self.dE_Ni
        elif foil == 'Cu':
            n = 10
            mol_mass_Cu = 63.546 #g/mol
            mass_density = np.array((22.338, 22.325, 22.313, 22.284, 22.443, 22.396, 22.320, 22.401, 22.425, 22.314)) /1e3
            sigma_mass_density = np.array((0.048, 0.028, 0.043, 0.027, 0.028, 0.012, 0.014, 0.033, 0.041, 0.047))/1e3
            mass_density = (mass_density*N_A)/mol_mass_Cu #nuclei/cm^2
            sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Cu
            E = self.E_Cu; dE = self.dE_Cu
            #print(mass_density[-1])

        return mass_density, sigma_mass_density, E, dE



    def cross_section(self,react_func, foil, filename, n, reaction, plot_CS=False):
        ###  example :  cross_section(Cu_57Ni, 'Cu', 'Cu_57Ni.csv/', 10)
        lamb = react_func[-1]
        mass_density, sigma_mass_density, E, dE = self.mass_density(foil)
        A0 = np.zeros(n)
        sigma_A0 = np.zeros(n)
        CS = np.zeros(n)
        sigma_CS = np.zeros(n)
        filename = self.path + filename


        #lamb, mass_density, sigma_mass_density = self.get_variables()

        for i in range(n):
            A0_val = np.genfromtxt(filename, delimiter=',', usecols=[i])
            A0[i] = A0_val[0]
            sigma_A0[i] = A0_val[1]


        for j in range(n):

            CS[j] = A0[j] / (mass_density[j] * self.I[j]*(1-np.exp(-lamb*self.irr_time)))*1e21   #mb  ###1 barn = 1e-24 cm^2


        if plot_CS:
            self.plot_CrossSections(E, CS, reaction)

        return A0, E, CS, self.I

    def plot_CrossSections(self, E, CS, reaction):
        plt.plot(E, CS, '.')
        #plt.plot(E, CS, '--', linewidth=0.2)
        plt.xlabel('Beam energy (MeV)')
        plt.ylabel('Cross section (mb)')
        plt.title('Cross section for reaction {}'.format(reaction))
        path = os.getcwd()
        plt.savefig(path + '/CrossSections/{}_CS.png'.format(reaction),dpi=300)
        plt.show()

    #def print_statements(self):








CS = CrossSections(files[25])
#path = os.getcwd() + '/activity_csv/'
#file = path + 'Ni_61Cu.csv'

#file = path + 'Cu_57Ni.csv/'
#x = CS.cross_section(Cu_57Ni(), 'Cu', 'Cu_57Ni.csv', 10)

x = CS.cross_section(Ir_193mPt(), 'Ir', 'Ir_193mPt.csv', 10, 'Ir_193mPt', plot_CS=True)
#x = CS.cross_section(Cu_57Ni(), 'Cu', 'Cu_57Ni.csv', 10, 'Cu_57Ni', plot_CS=True)
#x = CS.cross_section(Cu_64Cu(), 'Cu', 'Cu_64Cu.csv', 10, 'Cu_64Cu', plot_CS=True)
#print(x)

#CS.Get_CS(file, 10)
#CS.Cs(Ir_193mPt(), Ir_foil())
#CS.cross_section(file, 10)








"""
#from beam_current import *

path = os.getcwd() + '/activity_csv/'


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

def weighted_beam_current(E):
    return 1.54991205 * E + 90.99143333

def cross_section(filename):
    mass_density = Ni()
    n = 10
    A0_array = np.zeros(n)
    sigma = np.zeros(n)
    for i in range(n):
        A0 = np.genfromtxt(filename, delimiter=',', usecols=[i])
        #A0_array[i] = A0[0]
    #print(A0_array)
    #for j in range(n):
    #    sigma[j] = A0_array[j] / (mass_density[j] * I_weighted*(1-np.exp(-lamb*irr_time)))   #I_weighted is suppose to be weighed chi squared for beam current.
    #return sigma

#cross_section('Ni_57Co.csv')
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



CS = cross_section(path+'Ir_193mPt.csv')
#plot_CrossSections(E,CS, 'Ni_57Co')
"""
