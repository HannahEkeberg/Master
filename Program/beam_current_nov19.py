import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate
from scipy.constants import N_A, elementary_charge
import sys
from scipy.stats import norm
from scipy.optimize import curve_fit, minimize_scalar



from beam_current_FoilReact import *    #Program where info regarding foil reaction.
from ziegler_sorting import *  #sorting of ziegler list etc

"""
Ziegler energies and fluxes are used to first find energy and fluxes in different reactions
Function get_FWHM() returns FWHM for specific type of foils in a list for each foil.
Function beam_current() takes in foil and reaction, and returns the dI and I beamcurrent.
Uses functions data and E_flux_integral to get spline and integrals, from equation used
Function WABE (weighted average beam energy) returns the average energy for each foil, and dE
Finally plotting E,I with errorbars.

"""

from scipy.signal import chirp, find_peaks, peak_widths


#ziegler_file = '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/E_foils_fluxes.csv'
#filenames = [ziegler_file_SS_n10, ziegler_file_SS_n5,ziegler_file_SS_0, ziegler_file_SS_p5,ziegler_file_SS_p10, ziegler_file_Ni_n10]
#names = ['-SS-10%','-SS-5%', '-SS0%', '-SS+5%', '-SS+10%', '-Ni-10%']


#E_Ni, F_Ni, E_Cu, F_Cu, E_Fe, F_Fe = sort_ziegler(ziegler_file)
#name = 'Before variance minimization'#'-Ni+10%'


from variance_minimization import ziegler_files

files,names = ziegler_files()
n=0
name = names[n]; file = files[n]
E_Ni, F_Ni, E_Cu, F_Cu, E_Fe, F_Fe = sort_ziegler(file)




def get_sigmaE(E,F, foil):
    path_to_folder = os.getcwd()
    colors = ['mediumpurple', 'cyan', 'palevioletred', 'darkorange', 'forestgreen', 'orchid', 'dodgerblue', 'lime', 'crimson', 'indianred']
    dEr = np.zeros(len(E))
    dEl = np.zeros(len(E))
    for i in range(len(E)):
        M_F = np.max(F[i])  #max Flux
        Min_F = np.min(F[i])
        #mean_E = np.mean(F[i])
        #print(M_F)
        HM_F = 0.5*M_F      #Half max Flux

        def lin_interp(x, y, i, half):
            return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

        def half_max_F(E,F):
            half = max(F)/2.0
            signs = np.sign(np.add(F, -half))
            zero_crossings = (signs[0:-2] != signs[1:-1])
            zero_crossings_i = np.where(zero_crossings)[0]
            return [lin_interp(E, F, zero_crossings_i[0], half),
                    lin_interp(E, F, zero_crossings_i[1], half)]

        hmx = half_max_F(E[i], F[i])
        (mu,sigma) = norm.fit(E[i])
        fwhm = hmx[1]-hmx[0]

        #print(hmx[1], hmx[0])

        #print(mu-hmx[0], hmx[1]-mu)
        dEl[i] = mu-hmx[0]; dEr[i] = hmx[1]-mu
        half = max(F[i])/2.0
        #plt.plot(E[i], F[i], linewidth=0.6, color='navy')
        #plt.plot(hmx, [half, half], linewidth=0.8, color=colors[i], label='fwhm={0:.2f}'.format(fwhm))
        #plt.axhline(HM_F)
        #print(M_F, Min_F)
        #plt.vlines(mu, ymin=0.0, ymax = M_F, linewidth=0.4, linestyle='--')#, label=r'$\mu=${}'.format(mu))
    #plt.title('Energy distribution for {}-foils'.format(foil))
    #plt.xlabel('Energy, MeV')
    #plt.ylabel(r'Relative deuteron flux, $d\phi/dE$')
    #plt.legend()
    #plt.savefig(path_to_folder + '/BeamCurrent/' +foil+'_flux.png', dpi=300)
    #plt.show()
    return dEl, dEr


dEl_Ni, dEr_Ni = get_sigmaE(E_Ni,F_Ni, 'Ni')
dEl_Cu, dEr_Cu = get_sigmaE(E_Cu, F_Cu, 'Cu')
dEl_Fe, dEr_Fe = get_sigmaE(E_Fe, F_Fe, 'Fe')

def data(filename):    #Given the set of data points (x[i], y[i]) determine a smooth spline approximation
    E_mon = np.loadtxt(filename, usecols=[0], skiprows=6)
    Cs = np.loadtxt(filename, usecols=[1], skiprows=6)
    sigma_Cs = np.loadtxt(filename, usecols=[2], skiprows=6)

    tck = interpolate.splrep(E_mon, Cs, s=0)
    sigma_tck = interpolate.splrep(E_mon, sigma_Cs, s=0)
    return E_mon, Cs, sigma_Cs, tck, sigma_tck


def E_flux_integral(E_mon, Cs, sigma_Cs, tck, sigma_tck, E, F):  #integral in eq. int CS*d(phi)/dE dE = Cs*F dE
    reaction_integral = []
    uncertainty_integral = []

    for i in range(len(E)):
        Cs_ = interpolate.splev(E[i], tck, der=0)*1e-27 #mb--> 1e-27 cm^2. #gives interpolated cross section
        sigma_Cs_ = interpolate.splev(E[i], sigma_tck, der=0) * 1e-27
        relative_sigma_Cs = sigma_Cs_/ Cs_


        int_uncertainty = np.trapz(F[i]*relative_sigma_Cs, E[i])/np.trapz(F[i],E[i])
        uncertainty_integral.append(int_uncertainty)
        int_reaction = np.trapz(F[i]*Cs_, E[i])/np.trapz(F[i],E[i])
        reaction_integral.append(int_reaction)

    return uncertainty_integral, reaction_integral  #integral is now average percent uncertainty in cross section


def beam_current(foil, react):

    irr_time = 3600; sigma_irr_time = 3 # seconds, sigma is personally estimated
    integral = []


    if foil == 'Fe':
        F = F_Fe
        E = E_Fe  #from ziegler

        IAEA_Cs, A0, sigma_A0, lambda_, mass_density, sigma_mass_density = Fe_foil(react)
        #print(mass_density[0])
        E_mon, Cs, sigma_Cs, tck, sigma_tck = data(IAEA_Cs) #from monitor foils


    elif foil == 'Ni':
        IAEA_Cs, A0, sigma_A0, lambda_, mass_density, sigma_mass_density = Ni_foil(react)
        #print(mass_density[0])
        E_mon, Cs, sigma_Cs, tck, sigma_tck = data(IAEA_Cs) #from monitor foils
        F = F_Ni; E = E_Ni  #from ziegler

    elif foil == 'Cu':
        IAEA_Cs, A0, sigma_A0, lambda_, mass_density, sigma_mass_density = Cu_foil(react)
        E_mon, Cs, sigma_Cs, tck, sigma_tck = data(IAEA_Cs) #from monitor foils
        F = F_Cu; E = E_Cu  #from ziegler


    uncertainty_integral, reaction_integral = E_flux_integral(E_mon, Cs, sigma_Cs, tck, sigma_tck, E, F)
    uncertainty_integral = np.array((uncertainty_integral))

    I = A0 *elementary_charge*1e9 / (mass_density*(1-np.exp(-lambda_*irr_time))*reaction_integral)
    dI = I * np.sqrt((sigma_A0/A0)**2 + (sigma_mass_density/mass_density)**2 + (sigma_irr_time/irr_time)**2 + uncertainty_integral**2)

    I = np.array(I)
    return I, dI


def WABE(target): # Weighted Average Beam Energy

    if target == 'Ni':
        E = E_Ni; F = F_Ni
        dEl, dEr = dEl_Ni, dEr_Ni
    elif target == 'Cu':
        E = E_Cu; F = F_Cu
        dEl, dEr = dEl_Cu, dEr_Cu
    elif target == 'Fe':
        E = E_Fe; F = F_Fe
        dEl, dEr = dEl_Fe, dEr_Fe

    energy = []
    for index, item in enumerate(E):
        E[index] = np.array(E[index])
        E_int = np.trapz(F[index]*E[index], E[index])/np.trapz(F[index],E[index])
        energy.append(E_int)

    energy = np.array(energy)
    return energy, [dEl, dEr]

WE_Fe, sigma_WE_Fe = WABE('Fe')
WE_Ni, sigma_WE_Ni = WABE('Ni')
WE_Cu, sigma_WE_Cu = WABE('Cu')


I_Fe_56Co, dI_Fe_56Co = beam_current('Fe', 'Fe_56Co')
I_Ni_61Cu, dI_Ni_61Cu = beam_current('Ni', 'Ni_61Cu')
I_Ni_56Co, dI_Ni_56Co = beam_current('Ni', 'Ni_56Co')
I_Ni_58Co, dI_Ni_58Co = beam_current('Ni', 'Ni_58Co')
I_Cu_62Zn, dI_Cu_62Zn = beam_current('Cu', 'Cu_62Zn')    ##causes some problem in dI    RuntimeWarning: invalid value encountered in true_divide (in dI)
I_Cu_63Zn, dI_Cu_63Zn = beam_current('Cu', 'Cu_63Zn')    ##causes some problem in dI. Caused by I=0
I_Cu_65Zn, dI_Cu_65Zn = beam_current('Cu', 'Cu_65Zn')



def least_sqaures(data, model):
    return np.sum( (data-model)**2 / data)

def MSE(data, model):
    n = len(data)
    return 1/n * np.sum((data-model)**2)


def linear_fit():
    from sklearn import linear_model

    I = np.concatenate((I_Fe_56Co, I_Ni_61Cu, I_Ni_56Co, I_Ni_58Co, I_Cu_62Zn, I_Cu_63Zn, I_Cu_65Zn), axis=0).reshape(-1,1)
    E = np.concatenate((WE_Fe, WE_Ni, WE_Ni, WE_Ni, WE_Cu, WE_Cu, WE_Cu),axis=0).reshape(-1,1)
    #print(type(I), I.shape)
    I_new = I.T[0]
    E_new = E.T[0]

    index = np.where(I_new > 0.5)
    I_new = I_new[index].reshape(-1,1)
    E_new = E_new[index].reshape(-1,1)
    ##plt.plot(E_new, I_new, '.')
    #    plt.show()
    #I_new = I_new.reshape(-1,1)
    #I_new = I[I[i]!=0]


    #print(I_new)


    linreg = linear_model.LinearRegression().fit(E_new,I_new)
    I_pred = linreg.predict(E_new)

    beta = linreg.coef_
    intercept = linreg.intercept_
    #print('I(E) = ',beta,'*E + ',intercept)
    MSE_val = MSE(I_new, I_pred)
    #least_sqaures_val = least_sqaures(I, I_pred)
    #print(MSE_val)
    #print("MSE with {}:".format(file, MSE_val))
    #print("Least squares with {}:")
    print(MSE_val)
    #print(least_sqaures_val)
    return E_new, I_pred, MSE_val


#linear_fit()

def zero_to_nan(values):   #for those values which have zero activity in foil, make NaN and not plot. Here 62,63Zn
    #Replace every 0 with 'nan' and return a copy.
    x = [float('nan') if x==0 else x for x in values]
    index = ~(np.isnan(x))
    return x, index


def plot(name):

    E, I_pred = linear_fit()[:-1]
    plt.plot(E, I_pred, linewidth=1.0, color='red', label='Fit, mse = {0:.2f}'.format(linear_fit()[-1]))

    plt.axhline(128.5, linestyle='-.', linewidth=0.4, label='Monitor current 128.5 nA')


    plt.plot(WE_Fe,I_Fe_56Co, '.', label=r'$^{nat}$Fe(d,x)$^{56}$Co')
    plt.errorbar(WE_Fe, I_Fe_56Co, color='green', linewidth=0.001, xerr=sigma_WE_Fe, yerr=dI_Fe_56Co, elinewidth=0.5, ecolor='k', capthick=0.5 )

    plt.plot(WE_Ni,I_Ni_61Cu, '.', label=r'$^{nat}$Ni(d,x)$^{61}$Cu')
    plt.errorbar(WE_Ni, I_Ni_61Cu, color='green', linewidth=0.001, xerr=sigma_WE_Ni, yerr=dI_Ni_61Cu, elinewidth=0.5, ecolor='k', capthick=0.5 )

    plt.plot(WE_Ni,I_Ni_56Co, '.', label=r'$^{nat}$Ni(d,x)$^{56}$Co')
    plt.errorbar(WE_Ni, I_Ni_56Co, color='green', linewidth=0.001, xerr=sigma_WE_Ni, yerr=dI_Ni_56Co, elinewidth=0.5, ecolor='k', capthick=0.5 )

    plt.plot(WE_Ni,I_Ni_58Co, '.', label=r'$^{nat}$Ni(d,x)$^{58}$Co')
    plt.errorbar(WE_Ni, I_Ni_58Co, color='green', linewidth=0.001, xerr=sigma_WE_Ni, yerr=dI_Ni_58Co, elinewidth=0.5, ecolor='k', capthick=0.5 )

    #index = zero_to_nan(I_Cu_62Zn)[-1]
    #plt.plot(WE_Cu[index],I_Cu_62Zn[index], '.', label=r'$^{nat}$Cu(d,x)$^{62}$Zn')
    #plt.errorbar(WE_Cu[index], I_Cu_62Zn[index], color='green', linewidth=0.001, xerr=sigma_WE_Cu[index], yerr=dI_Cu_62Zn[index], elinewidth=0.5, ecolor='k', capthick=0.5 )
    plt.plot(WE_Cu,I_Cu_62Zn, '.', label=r'$^{nat}$Cu(d,x)$^{62}$Zn')
    plt.errorbar(WE_Cu, I_Cu_62Zn, color='green', linewidth=0.001, xerr=sigma_WE_Cu, yerr=dI_Cu_62Zn, elinewidth=0.5, ecolor='k', capthick=0.5 )

    #index = zero_to_nan(I_Cu_63Zn)[-1]
    #plt.plot(WE_Cu[index],I_Cu_63Zn[index], '.', label=r'$^{nat}$Cu(d,x)$^{63}$Zn')
    #plt.errorbar(WE_Cu[index], I_Cu_63Zn[index], color='green', linewidth=0.001, xerr=sigma_WE_Cu[index], yerr=dI_Cu_63Zn[index], elinewidth=0.5, ecolor='k', capthick=0.5 )
    plt.plot(WE_Cu,I_Cu_63Zn, '.', label=r'$^{nat}$Cu(d,x)$^{63}$Zn')
    plt.errorbar(WE_Cu, I_Cu_63Zn, color='green', linewidth=0.001, xerr=sigma_WE_Cu, yerr=dI_Cu_63Zn, elinewidth=0.5, ecolor='k', capthick=0.5 )

    plt.plot(WE_Cu,I_Cu_65Zn, '.', label=r'$^{nat}$Cu(d,x)$^{65}$Zn')
    plt.errorbar(WE_Cu, I_Cu_65Zn, color='green', linewidth=0.001, xerr=sigma_WE_Cu, yerr=dI_Cu_65Zn, elinewidth=0.5, ecolor='k', capthick=0.5 )
    #namez = names[file]
    plt.title('Beam current monitor foils {}'.format(name))
    plt.xlabel('Energy, MeV')
    plt.ylabel('Measured deuteron beam current, nA')
    plt.legend(fontsize='x-small')
    path = os.getcwd()
    plt.savefig(path+'/BeamCurrent/' + name +'.png', dpi=300)
    plt.show()


    #plot(names[file])
plot(name)
