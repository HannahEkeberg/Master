#ziegler_file = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_fluxes.csv'

#ziegler_file = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_fluxes.csv'
import numpy as np

def sort_ziegler(ziegler_file):
    ziegler_foil = np.loadtxt(ziegler_file, dtype="str", delimiter=",", usecols=[0], skiprows=1)
    ziegler_E = np.genfromtxt(ziegler_file, delimiter=',', usecols=[1], skip_header=1)
    ziegler_flux = np.genfromtxt(ziegler_file, delimiter=',', usecols=[2], skip_header=1)
    #arrays are 2818 long

    ziegler_foil = np.loadtxt(ziegler_file, dtype="str", delimiter=",", usecols=[0], skiprows=1)
    ziegler_E = np.genfromtxt(ziegler_file, delimiter=',', usecols=[1], skip_header=1)
    ziegler_flux = np.genfromtxt(ziegler_file, delimiter=',', usecols=[2], skip_header=1)

    E_Ni01=[]; E_Ni02=[]; E_Ni03=[]; E_Ni04=[]; E_Ni05=[]; E_Ni06=[]; E_Ni07=[]; E_Ni08=[]; E_Ni09=[]; E_Ni10=[]
    F_Ni01=[]; F_Ni02=[]; F_Ni03=[]; F_Ni04=[]; F_Ni05=[]; F_Ni06=[]; F_Ni07=[]; F_Ni08=[]; F_Ni09=[]; F_Ni10=[]

    E_Cu01=[]; E_Cu02=[]; E_Cu03=[]; E_Cu04=[]; E_Cu05=[]; E_Cu06=[]; E_Cu07=[]; E_Cu08=[]; E_Cu09=[]; E_Cu10=[]
    F_Cu01=[]; F_Cu02=[]; F_Cu03=[]; F_Cu04=[]; F_Cu05=[]; F_Cu06=[]; F_Cu07=[]; F_Cu08=[]; F_Cu09=[]; F_Cu10=[]

    E_Fe01=[]; E_Fe02=[]; E_Fe03=[]; E_Fe04=[]; E_Fe05=[]; E_Fe06=[]; E_Fe07=[]; E_Fe08=[]; E_Fe09=[]; E_Fe10=[]
    F_Fe01=[]; F_Fe02=[]; F_Fe03=[]; F_Fe04=[]; F_Fe05=[]; F_Fe06=[]; F_Fe07=[]; F_Fe08=[]; F_Fe09=[]; F_Fe10=[]



    for index, item in enumerate(ziegler_foil):
        if item == 'Ni01':
            E_Ni01.append(ziegler_E[index])
            F_Ni01.append(ziegler_flux[index])
        elif item == 'Ni02':
            E_Ni02.append(ziegler_E[index])
            F_Ni02.append(ziegler_flux[index])
        elif item == 'Ni03':
            E_Ni03.append(ziegler_E[index])
            F_Ni03.append(ziegler_flux[index])
        elif item == 'Ni04':
            E_Ni04.append(ziegler_E[index])
            F_Ni04.append(ziegler_flux[index])
        elif item == 'Ni05':
            E_Ni05.append(ziegler_E[index])
            F_Ni05.append(ziegler_flux[index])
        elif item == 'Ni06':
            E_Ni06.append(ziegler_E[index])
            F_Ni06.append(ziegler_flux[index])
        elif item == 'Ni07':
            E_Ni07.append(ziegler_E[index])
            F_Ni07.append(ziegler_flux[index])
        elif item == 'Ni08':
            E_Ni08.append(ziegler_E[index])
            F_Ni08.append(ziegler_flux[index])
        elif item == 'Ni09':
            E_Ni09.append(ziegler_E[index])
            F_Ni09.append(ziegler_flux[index])
        elif item == 'Ni10':
            E_Ni10.append(ziegler_E[index])
            F_Ni10.append(ziegler_flux[index])

        if item == 'Cu01':
            E_Cu01.append(ziegler_E[index])
            F_Cu01.append(ziegler_flux[index])
        elif item == 'Cu02':
            E_Cu02.append(ziegler_E[index])
            F_Cu02.append(ziegler_flux[index])
        elif item == 'Cu03':
            E_Cu03.append(ziegler_E[index])
            F_Cu03.append(ziegler_flux[index])
        elif item == 'Cu04':
            E_Cu04.append(ziegler_E[index])
            F_Cu04.append(ziegler_flux[index])
        elif item == 'Cu05':
            E_Cu05.append(ziegler_E[index])
            F_Cu05.append(ziegler_flux[index])
        elif item == 'Cu06':
            E_Cu06.append(ziegler_E[index])
            F_Cu06.append(ziegler_flux[index])
        elif item == 'Cu07':
            E_Cu07.append(ziegler_E[index])
            F_Cu07.append(ziegler_flux[index])
        elif item == 'Cu08':
            E_Cu08.append(ziegler_E[index])
            F_Cu08.append(ziegler_flux[index])
        elif item == 'Cu09':
            E_Cu09.append(ziegler_E[index])
            F_Cu09.append(ziegler_flux[index])
        elif item == 'Cu10':
            E_Cu10.append(ziegler_E[index])
            F_Cu10.append(ziegler_flux[index])

        if item == 'Fe01':
            E_Fe01.append(ziegler_E[index])
            F_Fe01.append(ziegler_flux[index])
        elif item == 'Fe02':
            E_Fe02.append(ziegler_E[index])
            F_Fe02.append(ziegler_flux[index])
        elif item == 'Fe03':
            E_Fe03.append(ziegler_E[index])
            F_Fe03.append(ziegler_flux[index])

    E_Ni = [E_Ni01, E_Ni02, E_Ni03, E_Ni04, E_Ni05, E_Ni06, E_Ni07, E_Ni08, E_Ni09, E_Ni10]
    F_Ni = [F_Ni01, F_Ni02, F_Ni03, F_Ni04, F_Ni05, F_Ni06, F_Ni07, F_Ni08, F_Ni09, F_Ni10]

    E_Cu = [E_Cu01, E_Cu02, E_Cu03, E_Cu04, E_Cu05, E_Cu06, E_Cu07, E_Cu08, E_Cu09, E_Cu10]
    F_Cu = [F_Cu01, F_Cu02, F_Cu03, F_Cu04, F_Cu05, F_Cu06, F_Cu07, F_Cu08, F_Cu09, F_Cu10]

    E_Fe = [E_Fe01, E_Fe02, E_Fe03]
    F_Fe = [F_Fe01, F_Fe02, F_Fe03]

    return E_Ni, F_Ni, E_Cu, F_Cu, E_Fe, F_Fe
