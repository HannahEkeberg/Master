import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate
from scipy.constants import N_A, elementary_charge
import sys
from scipy.stats import norm
from scipy.optimize import curve_fit, minimize_scalar

#from single_decay_A0 import *
from foil_info import *
import os

"""
Fe_56Co = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/_data_files/fed56cot/fed56cot.txt'
Ni_61Cu = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/nid61cut/nid61cut.txt'
Ni_56Co = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/nid56cot/nid56cot.txt'
Ni_58Co = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/nid58cot/nid58cot.txt'
Cu_62Zn = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/cud62znt/cud62znt.txt'
Cu_63Zn = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/cud63znt/cud63znt.txt'
Cu_65Zn = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/cud65znt/cud65znt.txt'
"""
path = os.getcwd() + '/activity_csv/'

def getA0(filename, n):# n=numb of foils. Get activity from csv to either foil function
    A0_array = np.zeros(n)
    sigma_A0_array = np.zeros(n)
    for i in range(n):
        A0 = np.genfromtxt(filename, delimiter=',', usecols=[i])
        A0_array[i] = A0[0]
        sigma_A0_array[i] = A0[1]
    return(A0_array, sigma_A0_array)

def Ir_foil():
    n=10
    mol_mass_Ir = 192.217 #g/mol
    mass_density = np.array((55.174, 55.601, 55.643, 56.000, 55.161, 55.731, 56.685, 58.030, 56.669, 55.065))/1e3
    sigma_mass_density = np.array((0.053, 0.238, 0.121, 0.109, 0.081, 0.088, 0.085, 0.130, 0.043, 0.055))/1e3
    mass_density = (mass_density*N_A)/mol_mass_Ir #nuclei/cm^2
    sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Ir

    return mass_density, sigma_mass_density


def Fe_foil(react): # returns A0, sigma_A0, lambda
    n = 3 #number of foils
    mol_mass_Fe = 55.8450 #g/mol
    mass_density = np.array((20.030, 20.017, 19.948))/1e3
    sigma_mass_density = np.array((0.110, 0.034, 0.114))/1e3
    mass_density = (mass_density*N_A)/mol_mass_Fe #nuclei/cm^2
    sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Fe

    if react == 'Fe_56Co':
        A_file = path + 'Fe_56Co.csv'
        A0, sigma_A0 = getA0(A_file, n)
        lambda_ = Fe_56Co()[-1]   #from foil_info.py
        IAEA_Cs = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/fed56cot/fed56cot.txt' ## IAEA mon data

    return IAEA_Cs, A0, sigma_A0, lambda_, mass_density, sigma_mass_density

def Ni_foil(react):
    n = 10
    mol_mass_Ni = 58.69340 #g/mol
    mass_density = np.array((22.772, 23.118,22.338,20.704,21.768,22.861,23.092,22.409,21.741,23.093))/1e3   #g/cm^2, divide by 1e3 to get g instead of mg.
    sigma_mass_density = np.array((0.138, 0.096, 0.066, 0.068, 0.045, 0.123,0.078, 0.124, 0.073, 0.024))/1e3 #g/cm^2
    mass_density = (mass_density*N_A)/mol_mass_Ni #nuclei/cm^2
    sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Ni

    if react == 'Ni_61Cu':
        A_file = path + 'Ni_61Cu.csv'
        A0, sigma_A0 = getA0(A_file, n)
        lambda_ = Ni_61Cu()[-1]   #from foil_info.py
        IAEA_Cs = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/nid61cut/nid61cut.txt'
        return IAEA_Cs, A0, sigma_A0, lambda_, mass_density, sigma_mass_density

    elif react == 'Ni_56Co':
        A_file = path + 'Ni_56Co.csv'
        A_56Ni_file = path + 'Ni_56Ni.csv'
        A0_56Co, sigma_A0_56Co = getA0(A_file, n)
        A0_56Ni, sigma_A0_56Ni = getA0(A_56Ni_file, n)


        lambda_56Co = Ni_56Co()[-1]   #from foil_info.py
        lambda_56Ni = Ni_56Ni()[-1]
        IAEA_Cs = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/nid56cot/nid56cot.txt'

        return IAEA_Cs, A0_56Co, sigma_A0_56Co, A0_56Ni, sigma_A0_56Ni, lambda_56Co, lambda_56Ni, mass_density, sigma_mass_density

    elif react == 'Ni_58Co':
        ###Cumulative cross section from IAEA
        A_file = path + 'Ni_58Co.csv'
        A_58mCo_file = path + 'Ni_58mCo.csv'
        A0_58Co, sigma_A0_58Co = getA0(A_file, n)
        A0_58mCo, sigma_A0_58mCo = getA0(A_58mCo_file, n)

        lambda_58Co = Ni_58Co()[-1]   #from foil_info.py
        lambda_58mCo = Ni_58Co()[1]   # Isomer decay with no gammas, 58mCo activities are calculated in Ni_58Co
        IAEA_Cs = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/nid58cot/nid58cot.txt'
        return IAEA_Cs, A0_58Co, sigma_A0_58Co, A0_58mCo, sigma_A0_58mCo, lambda_58Co, lambda_58mCo, mass_density, sigma_mass_density



def Cu_foil(react):
    n = 10
    mol_mass_Cu = 63.546 #g/mol
    mass_density = np.array((22.338, 22.325, 22.313, 22.284, 22.443, 22.396, 22.320, 22.401, 22.425, 22.314)) /1e3
    sigma_mass_density = np.array((0.048, 0.028, 0.043, 0.027, 0.028, 0.012, 0.014, 0.033, 0.041, 0.047))/1e3
    mass_density = (mass_density*N_A)/mol_mass_Cu #nuclei/cm^2
    sigma_mass_density = (sigma_mass_density*N_A)/mol_mass_Cu

    if react   == 'Cu_62Zn':
        A_file = path + 'Cu_62Zn.csv'
        A0, sigma_A0 = getA0(A_file, n)
        lambda_ = Cu_62Zn()[-1]   #from foil_info.py
        IAEA_Cs = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/cud62znt/cud62znt.txt'

    elif react == 'Cu_63Zn':
        A_file = path + 'Cu_63Zn.csv'
        A0, sigma_A0 = getA0(A_file, n)
        lambda_ = Cu_63Zn()[-1]   #from foil_info.py
        IAEA_Cs = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/cud63znt/cud63znt.txt'


    elif react == 'Cu_65Zn':
        A_file = path + 'Cu_65Zn.csv'
        A0, sigma_A0 = getA0(A_file, n)
        lambda_ = Cu_65Zn()[-1]   #from foil_info.py
        IAEA_Cs = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/monitor_data_files/cud65znt/cud65znt.txt'

    return IAEA_Cs, A0, sigma_A0, lambda_, mass_density, sigma_mass_density
