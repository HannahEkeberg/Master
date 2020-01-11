import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate
from scipy.constants import N_A, elementary_charge
import sys
from scipy.stats import norm
from scipy.optimize import curve_fit, minimize_scalar



from beam_current_FoilReact import *    #Program where info regarding foil reaction.
from ziegler_sorting import *  #sorting of ziegler list etc

from ZieglerFiles import ziegler_files

#files,names = ziegler_files()
#print(files)

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


class BeamCurrent:
    def __init__(self, ziegler_file, sort_ziegler, Fe_foil, Ni_Foil, Cu_foil):
        self.file = ziegler_file
        self.sort = sort_ziegler    # from ziegler_sorting.py
        self.Fe_foil = Fe_foil      # from beam_current_FoilReact
        self.Ni_foil = Ni_foil      # from beam_current_FoilReact
        self.Cu_foil = Cu_foil      # from beam_current_FoilReact
        self.E_Ni, self.F_Ni, self.E_Cu, self.F_Cu, self.E_Fe, self.F_Fe = self.sort(self.file)
        self.path = os.getcwd()

    def get_sigmaE(self, E, F, foil, makePlot=False):



        ### for cupper, flux is going up in end, delete those points!

        dEr = np.zeros(len(E))    # right uncertainty
        dEl = np.zeros(len(E))    # left uncertainty
        fwhm = np.zeros(len(E))
        half_max = []
        mu_array = np.zeros(len(E))
        for i in range(len(E)):
            M_F = np.max(F[i])     #max Flux
            Min_F = np.min(F[i])
            HM_F = 0.5*M_F         #Half max Flux

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
            half_max.append(hmx)
            (mu,sigma) = norm.fit(E[i])
            fwhm[i] = hmx[1]-hmx[0]
            dEl[i] = mu-hmx[0]; dEr[i] = hmx[1]-mu   #left and right uncertainty in energy
            mu_array[i] = mu
        if makePlot:
            self.Plot_energy_distribution(E,F,mu_array, fwhm, half_max, foil)  #make plot of energy distribution
        else:
            return dEl, dEr   #return left and right uncertainty

    def plot_distribution(self, foil, name):
        if foil == 'Ni':
            #print('way')
            self.get_sigmaE(self.E_Ni, self.F_Ni, foil, makePlot=True)
        elif foil == 'Cu':
            #print('way')
            self.get_sigmaE(self.E_Cu, self.F_Cu, foil, makePlot=True)
        elif foil == 'Fe':
            self.get_sigmaE(self.E_Fe, self.F_Fe, foil, makePlot=True)
            #print('way')
        path_to_folder = self.path + '/BeamCurrent/'
        plt.title('Energy distribution for {}-foils - '.format(foil) + name)
        plt.savefig(path_to_folder + foil + '_flux_distribution'+name+'.png', dpi=300)
        plt.show()


    def Plot_energy_distribution(self, E, F, mu, fwhm, half_max, foil):
        #print(half_max[0])
        #path_to_folder = self.path + '/BeamCurrent/'

        colors = ['mediumpurple', 'cyan', 'palevioletred', 'darkorange', 'forestgreen', 'orchid', 'dodgerblue', 'lime', 'crimson', 'indianred']
        for j in range(len(E)):
            plt.plot(E[j], F[j], color='navy', linewidth=0.7)
            half = np.max(F[j])*0.5
            plt.plot(half_max[j], [half, half], linewidth=0.8, color=colors[j], label='fwhm={0:.2f}'.format(fwhm[j]))
            plt.vlines(mu[j], ymin=0.0, ymax = np.max(F[j]), linewidth=0.4, linestyle='--')#, label=r'$\mu=${}'.format(mu))
        #plt.title('Energy distribution for {}-foils'.format(foil))
        plt.xlabel('Energy, MeV')
        plt.ylabel(r'Relative deuteron flux, $d\phi/dE$')
        plt.legend()
        #plt.savefig(path_to_folder + foil + '_flux_distribution'+name+'.png', dpi=300)
        #plt.show()


    def data(self, filename):   #Function that interpolates over energies and cross sections provided by IAEA
        E_mon = np.loadtxt(filename, usecols=[0], skiprows=6)
        Cs = np.loadtxt(filename, usecols=[1], skiprows=6)
        sigma_Cs = np.loadtxt(filename, usecols=[2], skiprows=6)

        tck = interpolate.splrep(E_mon, Cs, s=0)
        sigma_tck = interpolate.splrep(E_mon, sigma_Cs, s=0)
        return E_mon, Cs, sigma_Cs, tck, sigma_tck

    def E_flux_integral(self, Cs, sigma_Cs, tck, sigma_tck, E, F):
        reaction_integral = []    # from beam current equation
        uncertainty_integral = [] # make relative uncertainty, by interpolating over ziegler E
        for i in range(len(E)):
            Cs_ = interpolate.splev(E[i], tck, der=0)*1e-27 #mb--> 1e-27 cm^2. #gives interpolated cross section
            sigma_Cs_ = interpolate.splev(E[i], sigma_tck, der=0) * 1e-27
            relative_sigma_Cs = sigma_Cs_/ Cs_

            #print("relative sigma Cs:", relative_sigma_Cs)
            int_uncertainty = np.trapz(F[i]*relative_sigma_Cs, E[i])/np.trapz(F[i],E[i])
            uncertainty_integral.append(int_uncertainty)

            int_reaction = np.trapz(F[i]*Cs_, E[i])/np.trapz(F[i],E[i])
            reaction_integral.append(int_reaction)

        return uncertainty_integral, reaction_integral


    def calculate_beam_current(self, foil, react, print_terms=False):   # For non-cumulative cross sections from IAEA.
        irr_time = 3600; sigma_irr_time = 3 #seconds

        if foil == 'Fe':
            F = self.F_Fe
            E = self.E_Fe

            IAEA_Cs, A0, sigma_A0, lambda_, mass_density, sigma_mass_density = self.Fe_foil(react)  #from beam_current_FoilReact
            E_mon, Cs, sigma_Cs, tck, sigma_tck = self.data(IAEA_Cs) #from monitor foils


        elif foil == 'Ni':
            F = self.F_Ni
            E = self.E_Ni
            if react == 'Ni_61Cu':
                IAEA_Cs, A0, sigma_A0, lambda_, mass_density, sigma_mass_density = self.Ni_foil(react)  #from beam_current_FoilReact
            else:  # Cumulative activities
                IAEA_Cs, A0_dir, sigma_A0_dir, A0_nondir, sigma_A0_nondir, lambda_dir, lambda_nondir, mass_density, sigma_mass_density = self.Ni_foil(react)

            E_mon, Cs, sigma_Cs, tck, sigma_tck = self.data(IAEA_Cs) #from monitor foils

        elif foil == 'Cu':
            F = self.F_Cu
            E = self.E_Cu

            IAEA_Cs, A0, sigma_A0, lambda_, mass_density, sigma_mass_density = self.Cu_foil(react)  #from beam_current_FoilReact
            E_mon, Cs, sigma_Cs, tck, sigma_tck = self.data(IAEA_Cs) #from monitor foils


        uncertainty_integral, reaction_integral = self.E_flux_integral(Cs, sigma_Cs, tck, sigma_tck, E, F)
        uncertainty_integral = np.array((uncertainty_integral))

        if react=='Ni_58Co' or react=='Ni_56Co':
            BR = 1.0
            sigma_BR = 0.0 # no BR for these reactions.

            I = 1/(mass_density * reaction_integral)  * ( A0_dir*elementary_charge*1e9/(1-np.exp(-lambda_dir*irr_time)) + BR*A0_nondir*elementary_charge*1e9/ (1-np.exp(-lambda_nondir*irr_time)) )
            dI = I * np.sqrt((sigma_A0_dir/A0_dir)**2 + (sigma_A0_nondir/A0_nondir)**2  + (sigma_BR/BR)**2 + (sigma_mass_density/mass_density)**2 + (sigma_irr_time/irr_time)**2 + uncertainty_integral**2)

            if print_terms:
                print("mass density", mass_density)
                print("reaction integral", reaction_integral)
                print("underctainty integral", uncertainty_integral)
                print("activity direct", A0_dir)
                print("activity non-direct", A0_nondir)

        else:
            I = A0 *elementary_charge*1e9 / (mass_density*(1-np.exp(-lambda_*irr_time))*reaction_integral)
            dI = I * np.sqrt((sigma_A0/A0)**2 + (sigma_mass_density/mass_density)**2 + (sigma_irr_time/irr_time)**2 + uncertainty_integral**2)

        if print_terms:
            print("mass density", mass_density)
            print("reaction integral", reaction_integral)
            print("underctainty integral", uncertainty_integral)
            print("activity", A0)

        I = np.array(I)
        #print(react)
        #print(I)
        return I, dI


    def WABE(self, foil): #weighted average beam energy

        if foil == 'Ni':
            E = self.E_Ni; F  = self.F_Ni
        elif foil == 'Cu':
            E = self.E_Cu; F  = self.F_Cu
        elif foil == 'Fe':
            E = self.E_Fe; F  = self.F_Fe

        dEl, dEr = self.get_sigmaE(E, F, foil)
        energy = []
        for index, item in enumerate(E):
            E[index] = np.array(E[index])
            E_int = np.trapz(F[index]*E[index], E[index])/np.trapz(F[index],E[index])
            energy.append(E_int)

        energy=np.array(energy)

        return energy, [dEl, dEr]

    def specified_currents(self, uncertainty=False):
        I_Fe_56Co, dI_Fe_56Co = self.calculate_beam_current('Fe', 'Fe_56Co')
        I_Ni_61Cu, dI_Ni_61Cu = self.calculate_beam_current('Ni', 'Ni_61Cu')
        I_Ni_56Co, dI_Ni_56Co = self.calculate_beam_current('Ni', 'Ni_56Co')
        I_Ni_58Co, dI_Ni_58Co = self.calculate_beam_current('Ni', 'Ni_58Co')
        I_Cu_62Zn, dI_Cu_62Zn = self.calculate_beam_current('Cu', 'Cu_62Zn')    ##causes some problem in dI    RuntimeWarning: invalid value encountered in true_divide (in dI)
        I_Cu_63Zn, dI_Cu_63Zn = self.calculate_beam_current('Cu', 'Cu_63Zn')    ##causes some problem in dI. Caused by I=0
        I_Cu_65Zn, dI_Cu_65Zn = self.calculate_beam_current('Cu', 'Cu_65Zn')
        if uncertainty:
            return dI_Fe_56Co, dI_Ni_61Cu, dI_Ni_56Co, dI_Ni_58Co, dI_Cu_62Zn, dI_Cu_63Zn, dI_Cu_65Zn
        else:
            return I_Fe_56Co, I_Ni_61Cu, I_Ni_56Co, I_Ni_58Co, I_Cu_62Zn, I_Cu_63Zn, I_Cu_65Zn


    def specified_energies(self, uncertainty=False):
        E_Fe, dE_Fe = self.WABE('Fe')
        E_Ni, dE_Ni = self.WABE('Ni')
        E_Cu, dE_Cu = self.WABE('Cu')

        if uncertainty:
            return dE_Fe, dE_Ni, dE_Cu
        else:
            return E_Fe, E_Ni, E_Cu

        #print(I_Fe_56Co)


    def variance_minimization(self, compartment, name):
        # Compartment means foil number positions
        #variance minimization & standard deviation
        I_Fe_56Co, I_Ni_61Cu, I_Ni_56Co, I_Ni_58Co, I_Cu_62Zn, I_Cu_63Zn, I_Cu_65Zn = self.specified_currents()
        dI_Fe_56Co, dI_Ni_61Cu, dI_Ni_56Co, dI_Ni_58Co, dI_Cu_62Zn, dI_Cu_63Zn, dI_Cu_65Zn = self.specified_currents(uncertainty=True)
        WE_Fe, WE_Ni, WE_Cu = self.specified_energies()
        #print(WE_Fe, WE_Ni, WE_Cu)
        sigma_WE_Fe, sigma_WE_Ni, sigma_WE_Cu = self.specified_energies(uncertainty=True)

        I_Ni_61Cu = I_Ni_61Cu[compartment]; dI_Ni_61Cu = dI_Ni_61Cu[compartment]
        I_Ni_56Co = I_Ni_56Co[compartment]; dI_Ni_56Co = dI_Ni_56Co[compartment]
        I_Ni_58Co = I_Ni_58Co[compartment]; dI_Ni_58Co = dI_Ni_58Co[compartment]
        I_Cu_62Zn = I_Cu_62Zn[compartment]; dI_Cu_62Zn = dI_Cu_62Zn[compartment]
        I_Cu_63Zn = I_Cu_63Zn[compartment]; dI_Cu_63Zn = dI_Cu_63Zn[compartment]
        I_Cu_65Zn = I_Cu_65Zn[compartment]; dI_Cu_65Zn = dI_Cu_65Zn[compartment]


        WE_Ni = WE_Ni[compartment]
        WE_Cu = WE_Cu[compartment]
        dWE_Ni = np.array(([sigma_WE_Ni[0][compartment]], [sigma_WE_Ni[1][compartment]]))
        dWE_Cu = np.array(([sigma_WE_Cu[0][compartment]], [sigma_WE_Cu[1][compartment]] ))

        I_Ni = np.array((I_Ni_61Cu, I_Ni_56Co, I_Ni_58Co))
        dI_Ni = np.array((dI_Ni_61Cu, dI_Ni_56Co, dI_Ni_58Co))
        I_Cu = np.array((I_Cu_62Zn, I_Cu_63Zn, I_Cu_65Zn))
        dI_Cu = np.array((dI_Cu_62Zn, dI_Cu_63Zn, dI_Cu_65Zn))
        E_Ni = np.ones(len(I_Ni)) * WE_Ni
        E_Cu = np.ones(len(I_Cu)) * WE_Cu
        I = np.concatenate((I_Ni, I_Cu))
        dI = np.concatenate((dI_Ni, dI_Cu))
        E = np.concatenate((E_Ni, E_Cu))

        def I_model(x,b):
            # Linear model, with slope set to 0.
            # Since nothing between foil Cu07 and Ni07, the current degradation=0.
            m = 0
            b = np.ones(len(x))*b
            return x*m + b

        if compartment <= len(I_Fe_56Co):
            I_Fe = I_Fe_56Co[compartment]; dI_Fe = dI_Fe_56Co[compartment]
            WE_Fe = WE_Fe[compartment]

            I = np.append(I, I_Fe)
            dI = np.append(dI, dI_Fe)
            E = np.append(E, WE_Fe)
            dWE_Fe = np.array(([sigma_WE_Fe[0][compartment]], [sigma_WE_Fe[1][compartment]] ))



        popt, pcov = curve_fit(I_model, E, I, p0=128, sigma=dI, absolute_sigma=True)
        sigma_I_est = float( np.sqrt(np.diagonal(pcov)) ) #Uncertainty in the fitting parameters
        I_est = popt[0]

        chi_sq = self.chi_sqaured(I, I_est, sigma_I_est)
        #print(chi_sq)
        """
        xplot = np.linspace(min(E)-0.5, max(E)+0.5,len(E))
        plt.plot(xplot, I_model(E, *popt), label=r'Linear fit I={:.2f} $\pm$ {:.2f} nA'.format(I_est, sigma_I_est), linestyle='--', color='red')
        plt.plot(xplot,I_model(E,*(popt+sigma_I_est)), color='blue', linewidth=0.4, linestyle='-.')
        plt.plot(xplot,I_model(E,*(popt-sigma_I_est)), color='blue', linewidth=0.4,linestyle='-.', label=r'Uncertainty in fit, 1$\sigma$')

        plt.plot(WE_Ni, I_Ni_61Cu, marker='o', label=r'$Ni(d,x)^{61}Cu$')
        plt.errorbar(WE_Ni, I_Ni_61Cu, color='green', linewidth=0.001,xerr=dWE_Ni, yerr=dI_Ni_61Cu, elinewidth=0.5, ecolor='k', capthick=0.5 )
        plt.plot(WE_Ni, I_Ni_56Co, marker='o', label=r'$Ni(d,x)^{56}Co$')
        plt.errorbar(WE_Ni, I_Ni_56Co, color='green', linewidth=0.001,xerr=dWE_Ni, yerr=dI_Ni_56Co, elinewidth=0.5, ecolor='k', capthick=0.5 )
        plt.plot(WE_Ni, I_Ni_58Co, marker='o', label=r'$Ni(d,x)^{58}Co$')
        plt.errorbar(WE_Ni, I_Ni_58Co, color='green', linewidth=0.001,xerr=dWE_Ni, yerr=dI_Ni_58Co, elinewidth=0.5, ecolor='k', capthick=0.5 )


        #plt.plot(WE_Ni, I_Ni_61Cu‚ '.', label='Ni')

        plt.plot(WE_Cu, I_Cu_62Zn, marker='o', label=r'$Cu(d,x)^{62}Zn$')
        plt.errorbar(WE_Cu, I_Cu_62Zn, color='green', linewidth=0.001,xerr=dWE_Cu, yerr=dI_Cu_62Zn, elinewidth=0.5, ecolor='k', capthick=0.5 )
        plt.plot(WE_Cu, I_Cu_63Zn, marker='o', label=r'$Cu(d,x)^{63}Zn$')
        plt.errorbar(WE_Cu, I_Cu_63Zn, color='green', linewidth=0.001,xerr=dWE_Cu, yerr=dI_Cu_63Zn, elinewidth=0.5, ecolor='k', capthick=0.5 )
        plt.plot(WE_Cu, I_Cu_65Zn, marker='o', label=r'$Cu(d,x)^{65}Zn$')
        plt.errorbar(WE_Cu, I_Cu_65Zn, color='green', linewidth=0.001,xerr=dWE_Cu, yerr=dI_Cu_65Zn, elinewidth=0.5, ecolor='k', capthick=0.5 )

        if compartment <= len(I_Fe_56Co):
            #print("Yo")
            #print(type(WE_Fe), type(I_Fe_56Co))
            plt.plot(WE_Fe, I_Fe,marker='o', label=r'$Fe(d,x)^{56}Co$')
            plt.errorbar(WE_Fe, I_Fe, color='green', linewidth=0.001,xerr=dWE_Fe, yerr=dI_Fe, elinewidth=0.5, ecolor='k', capthick=0.5 )
            #plt.plot(E_Fe, I_Fe, 'o', label='Fe')
            #plt.errorbar(E_Fe, I_Fe, color='green', linewidth=0.001, yerr=dI_Fe, elinewidth=0.5, ecolor='k', capthick=0.5 )


        plt.xlabel('Energy, MeV')
        plt.ylabel('Beam Current, nA')
        plt.title(r'Linear fit for foils compartment {} - {}. $\chi^2$={:.2f} '.format(compartment+1, name, chi_sq))
        plt.legend(fontsize='x-small')
        plt.savefig('BeamCurrent/Chi_minimization/minimization_{}_'.format(compartment+1)+name+'.png', dpi=300)
        plt.show()
        """

        return WE_Ni, chi_sq



    def chi_sqaured(self, data, model, error): #error = stdv
        #dof = 5
        #return np.sum( (data - model/error )**2 )#( / (len(data-dof))
        n = len(data)
        #model = np.ones(n)*model
        #x = np.sum(data-model)
        #print(x)
        stdv = np.sqrt(1./n *  (np.sum(data-model)**2))
        #print(stdv)
        #print(stdv)
        return np.sum((data-model)**2 / stdv)
        #pass

    def linear_fit(self, makePlot=False):
        from sklearn import linear_model
        I = self.specified_currents()
        I = np.concatenate(I)
        WE = self.specified_energies()
        E = np.array((WE[0], WE[1], WE[1], WE[1], WE[-1], WE[-1], WE[-1]))
        E = np.concatenate(E)

        index = np.where (I > 0)
        I = I[index].reshape(-1,1)
        E = E[index].reshape(-1,1)

        linreg = linear_model.LinearRegression().fit(I,E)
        I_pred = linreg.predict(E)

        beta = linreg.coef_
        intercept = linreg.intercept_

        chi = self.least_squares(I, I_pred)


        if makePlot:
            plt.plot(E, I_pred, label='test')
            plt.show()



    def CurrentPlot(self, name):
        I_Fe_56Co, I_Ni_61Cu, I_Ni_56Co, I_Ni_58Co, I_Cu_62Zn, I_Cu_63Zn, I_Cu_65Zn = self.specified_currents()
        dI_Fe_56Co, dI_Ni_61Cu, dI_Ni_56Co, dI_Ni_58Co, dI_Cu_62Zn, dI_Cu_63Zn, dI_Cu_65Zn = self.specified_currents(uncertainty=True)
        WE_Fe, WE_Ni, WE_Cu = self.specified_energies()
        sigma_WE_Fe, sigma_WE_Ni, sigma_WE_Cu = self.specified_energies(uncertainty=True)



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
        plt.title('Beam current monitor foils - {}'.format(name))
        plt.xlabel('Energy, MeV')
        plt.ylabel('Measured deuteron beam current, nA')
        plt.ylim(0,350)
        plt.xlim(0,35)
        #ylim(top=350)
        #ylim(bottom=0)
        plt.legend(fontsize='x-small')
        path = os.getcwd()
        plt.savefig(path + '/BeamCurrent/beamcurrent_' + name  + '.png',dpi=300)
        #plt.savefig(path+'/BeamCurrent/' + name +'.png', dpi=300)
        plt.show()


files,names = ziegler_files()
ziegler_file = files[1]
#name = names[1]
#myclass = BeamCurrent(ziegler_file, sort_ziegler, Fe_foil, Ni_foil, Cu_foil)
#myclass.variance_minimization(7, 'SS_-5%')
#myclass.calculate_beam_current('Cu', 'Cu_65Zn', print_terms=True)
#myclass.CurrentPlot(name)

def run_all(compartment):
    chi_squared = []
    Ni_energies = []
    for i in range(len(files)):
        ### Remember that foil1 has index i=0.
        myclass = BeamCurrent(files[i], sort_ziegler, Fe_foil, Ni_foil, Cu_foil)
        #myclass.calculate_beam_current('Cu', 'Cu_65Zn', print_terms=True)
        #myclass.CurrentPlot(names[i])
        WE_Ni,chi_sq = myclass.variance_minimization(compartment, names[i])
        chi_squared.append(chi_sq)
        Ni_energies.append(WE_Ni)
    #print(Ni_energies, chi_squared)
    plt.plot(Ni_energies, chi_squared)
    plt.title(r'$\chi^2$')
    plt.xlabel('Deuteron energy entering {}th stack compartment (MeV)'.format(compartment+1))
    #plt.ylabel('')
    plt.savefig('BeamCurrent/Chi_minimization/chi_squared_comp_{}'.format(compartment+1))
    plt.show()




        #np.savetxt("{}.csv".format(save_results_to +  reaction_parent), np.array((A0_parent, sigma_A0_parent)), delimiter=",")

    #myclass.WABE()

    #myclass.plot_distribution('Ni', names[i])
    #myclass.plot_distribution('Cu', names[i])
    #myclass.plot_distribution('Fe', names[i])

run_all(8)
#path= '/Users/hannah/Documents/UIO/Masteroppgaven/Ziegler/'
#ziegler_file = path + 'E_foils_fluxes.csv'
#name = 'NoScalingParameter'
#myclass = BeamCurrent(ziegler_file, sort_ziegler, Fe_foil, Ni_foil, Cu_foil)
#myclass.plot_distribution('Cu', name)
#myclass.plot_distribution('Ni', name)
#myclass.plot_distribution('Fe', name)
#myclass.CurrentPlot(name)

#ziegler_file = path + 'E_foils_Fe_-1_fluxes.csv'
#name = 'Fe_-1'


#myclass.calculate_beam_current('Ni', 'Ni_61Cu')
#myclass.plot_distribution('Cu')
#myclass.calculate_beam_current('Cu', 'Cu_63Zn')
#myclass.WABE('Cu')
#myclass.CurrentPlot(name)
#myclass.linear_fit(makePlot=True)





"""

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
        dEl[i] = mu-hmx[0]; dEr[i] = hmx[1]-mu   #left and right uncertainty in energy
        #half = max(F[i])/2.0
        #plt.plot(E[i], F[i], color='navy', linewidth=0.7)
        #plt.plot(hmx, [half, half], linewidth=0.8, color=colors[i], label='fwhm={0:.2f}'.format(fwhm))
        #plt.axhline(HM_F)
        #print(M_F, Min_F)
        #        plt.vlines(mu, ymin=0.0, ymax = M_F, linewidth=0.4, linestyle='--')#, label=r'$\mu=${}'.format(mu))
    #plt.title('Energy distribution for {}-foils'.format(foil))
    #plt.xlabel('Energy, MeV')
    #plt.ylabel(r'Relative deuteron flux, $d\phi/dE$')
    #plt.legend()
    #plt.savefig(path_to_folder + '/BeamCurrent/' +foil+'_flux_distribution.png', dpi=300)
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
    print('I(E) = ',beta,'*E + ',intercept)
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
#plot(name)
"""
