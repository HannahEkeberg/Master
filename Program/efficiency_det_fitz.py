import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
#sys.path.insert(0, '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Calibration/Programmateriale/Detector_Info/')
#sys.path.insert(0,'/Users/hannah/Documents/UIO/Masteroppgaven/Data/Calibration/Class_eff_cal/')
#sys.path.insert(0, '/Users/hannah/⁨Dokumenter⁩/UIO⁩/Masteroppgaven⁩/⁨Data⁩/Calibration/Calibration_files/')


from detector_info_fitz import Detector_Information_fitz



class Efficiency_calculations(Detector_Information_fitz):

    def __init__(self, detector):
        Detector_Information_fitz.__init__(self) #Same as super, another option
        self.detector=detector

        #want to make a python file, which when imported the necessary information about the detector,
        #such as cal files (f_Cs137, f_Ba133, f_Eu152), live time, real time, and dates of spectra will be assigned automatically in class


        #getting information from detector_info.py
        list_detector = ['HPGE1_10','HPGE1_30','HPGE2_10','HPGE2_30', 'IDM1_53','IDM2_32','IDM3_40','IDM4_25', 'room131_5', 'room131_10', 'room131_15', 'room131_18', 'room131_22', 'room131_30', 'room131_40', 'room131_50', 'room131_60'] # ['HPGE', 'HPGE2', 'IDM1', 'IDM2', 'IDM3', 'IDM4', 'room_131_5cm', 'room_131_10cm', 'room_131_15cm', 'room_131_18cm', 'room_131_22cm', 'room_131_30cm', 'room_131_40cm', 'room_131_50cm', 'room_131_60cm']
        for element in list_detector:
            if element == self.detector:

                detector_function_from_class = getattr(Detector_Information_fitz, self.detector)
                self.live_time, self.time_delay, self.B, f_Cs137, f_Ba133, f_Eu152 = detector_function_from_class(self)

        #General information about each source
        self.t_half=np.array((30.08,10.551,13.517,2.6018))*365*24*3600 # seconds, Cs137, Ba133, Eu152, Na22
        self.t_half_uncert = np.array((0.09,0.011, 0.014,0.0022))*365*24*3600

        self.lamb=np.log(2)/self.t_half
        self.sigma_lamb= self.lamb*self.t_half_uncert/self.t_half

        self.A_0=np.array((38.55,39.89,39.29))*1e3 # Source activity Bq, Cs,Ba,Eu,Na
        self.sigma_A_0=0.3  #uncertainty in activity sources

        self.I_Cs137 = np.genfromtxt('Intensities_cal_sources.txt', usecols=[1], skip_header=3, skip_footer=30)
        self.I_Ba133 = np.genfromtxt('Intensities_cal_sources.txt', usecols=[1], skip_header=9, skip_footer=21)
        self.I_Eu152 = np.genfromtxt('Intensities_cal_sources.txt', usecols=[1], skip_header=19, skip_footer=0)

        self.sigma_I_Cs137 = np.genfromtxt('Intensities_cal_sources.txt', usecols=[2], skip_header=3, skip_footer=30)
        self.sigma_I_Ba133 = np.genfromtxt('Intensities_cal_sources.txt', usecols=[2], skip_header=9, skip_footer=21)
        self.sigma_I_Eu152 = np.genfromtxt('Intensities_cal_sources.txt', usecols=[2], skip_header=19, skip_footer=0)

        @staticmethod
        def information_from_files(filename):
            channel=np.loadtxt(filename, skiprows=37,usecols=[1])
            Nc=np.loadtxt(filename, usecols=[5], skiprows=37) #net area
            sigma_Nc=np.loadtxt(filename, usecols=[6], skiprows=37)
            E=np.loadtxt(filename, usecols=[0], skiprows=37)

            return E, Nc, sigma_Nc

        self.E_Cs, self.Nc_Cs, self.sigma_Nc_Cs = information_from_files.__func__(f_Cs137)
        self.E_Ba, self.Nc_Ba, self.sigma_Nc_Ba = information_from_files.__func__(f_Ba133)
        self.E_Eu, self.Nc_Eu, self.sigma_Nc_Eu = information_from_files.__func__(f_Eu152)
        #print(self.Nc_Cs)
        #for index, element in enumerate(self.Nc_Cs, start=0):
        #    if element == 0:
        #        del self.Nc_Cs, self.I_Cs[index], self.E_Cs[index], self.sigma_Nc_Cs[index]
        #print(self.Nc_Cs)

        self.eps_Cs = self.Nc_Cs*self.lamb[0]/(self.A_0[0]*self.I_Cs137*(1-np.exp(-self.lamb[0]*self.live_time[0]))*np.exp(-self.lamb[0]*self.time_delay[0]))
        self.eps_Ba = self.Nc_Ba*self.lamb[1]/(self.A_0[1]*self.I_Ba133*(1-np.exp(-self.lamb[1]*self.live_time[1]))*np.exp(-self.lamb[1]*self.time_delay[1]))
        self.eps_Eu = self.Nc_Eu*self.lamb[2]/(self.A_0[2]*self.I_Eu152*(1-np.exp(-self.lamb[2]*self.live_time[2]))*np.exp(-self.lamb[2]*self.time_delay[2]))

        print("Nc:{} ".format(self.Nc_Ba))
        print("lambda {}".format(self.lamb[1]))
        print("A0 {} ".format(self.A_0[1]))
        print("I {}".format(self.I_Ba133))
        print("live time {}, time delay {}".format(self.live_time, self.time_delay))

        self.sigma_eps_Cs = self.eps_Cs*(np.sqrt(self.sigma_lamb[0]/self.lamb[0])**2 + np.sqrt(self.sigma_Nc_Cs/self.Nc_Cs)**2 + np.sqrt(self.sigma_A_0/self.A_0[0])**2 + np.sqrt(self.sigma_I_Cs137/self.I_Cs137)**2)
        self.sigma_eps_Ba = self.eps_Ba*(np.sqrt(self.sigma_lamb[1]/self.lamb[1])**2 + np.sqrt(self.sigma_Nc_Ba/self.Nc_Ba)**2 + np.sqrt(self.sigma_A_0/self.A_0[1])**2 + np.sqrt(self.sigma_I_Ba133/self.I_Ba133)**2)
        self.sigma_eps_Eu = self.eps_Eu*(np.sqrt(self.sigma_lamb[2]/self.lamb[0])**2 + np.sqrt(self.sigma_Nc_Eu/self.Nc_Eu)**2 + np.sqrt(self.sigma_A_0/self.A_0[2])**2 + np.sqrt(self.sigma_I_Eu152/self.I_Eu152)**2)

        self.E = np.concatenate((self.E_Cs, self.E_Ba, self.E_Eu))
        self.eps = np.concatenate((self.eps_Cs, self.eps_Ba, self.eps_Eu))
        #self.sigma_E = np.concatenate((self.sigma_E_Cs, self.sigma_E_Ba, self.sigma_E_Eu))  not runnning with uncertainties in x
        self.sigma_eps = np.concatenate((self.sigma_eps_Cs, self.sigma_eps_Ba, self.sigma_eps_Eu))


    def efficiency_estimated_popt(self, E_gamma, B0, B1, B2, B3, B4):   #didnt work in curve fit so had to adapt. Function under is used in other subclasses.
        eps_est = B0*np.exp(-B1*E_gamma**B2)*(1-np.exp(-B3*E_gamma**B4))
        return eps_est

    def efficiency_estimated(self, E_gamma):
        eps_est = self.B[0]*np.exp(-self.B[1]*E_gamma**self.B[2])*(1-np.exp(-self.B[3]*E_gamma**self.B[4]))
        return eps_est



    def plot_data(self):
        plt.plot(self.E_Cs, self.eps_Cs, '.', color='magenta')
        plt.plot(self.E_Ba, self.eps_Ba,'.', color='green')
        plt.plot(self.E_Eu, self.eps_Eu, '.', color='blue')

        plt.errorbar(self.E_Cs, self.eps_Cs, color='green', linewidth=0.001, yerr=self.sigma_eps_Cs, elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        plt.errorbar(self.E_Ba, self.eps_Ba, color='magenta', linewidth=0.001, yerr=self.sigma_eps_Ba, elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        plt.errorbar(self.E_Eu, self.eps_Eu,color='blue', linewidth=0.001, yerr=self.sigma_eps_Eu, elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        plt.title("Datapoints and errors for calibration sources measured at {}".format(self.detector))
        plt.xlabel("Energy (keV)")
        plt.ylabel("Efficiency")
        plt.legend(["$^{137}Cs$", "$^{133}Ba$", "$^{152}Eu$"], loc='best')
        plt.show()

    def plot_func(self):

        B0=self.B[0]; B1=self.B[1];B2=self.B[2]; B3=self.B[3];B4=self.B[4]

        efficiency_estimated_popt = self.efficiency_estimated_popt(self.E,B0, B1,B2,B3,B4)

        #self.E = self.E[2:]; self.eps = self.eps[2:]; self.sigma_eps= self.sigma_eps[2:]
        index = ~(np.isnan(self.eps) | np.isnan(self.sigma_eps))  #if either eps OR sigma eps is NaN
        print(index)
        #print(self.eps)
        #print(self.sigma_eps)
        #print(self.eps[index])
        self.E = self.E[index]; self.eps = self.eps[index]; self.sigma_eps = self.sigma_eps[index]
        #self.eps_Cs = self.eps_Cs[index];self.eps_Ba = self.eps_Ba[index];self.eps_Eu = self.eps_Eu[index]
        #self.E_Cs = self.E_Cs[index];self.E_Ba = self.E_Ba[index];self.E_Eu = self.E_Eu[index]
        #self.sigma_eps_Cs = self.sigma_eps_Cs[index];self.sigma_eps_Ba = self.sigma_eps_Ba[index];self.sigma_eps_Eu = self.sigma_eps_Eu[index]

        print("E: {}".format(self.E))
        popt, pcov=curve_fit(self.efficiency_estimated_popt, self.E, self.eps, p0=np.array([B0, B1, B2, B3, B4]), sigma=self.sigma_eps, absolute_sigma=True, maxfev=1000000)#self.y_err/self.y) #John's numbers
        #popt = np.array((0.078, 2.0, 0.1637, 5.29e-5, 2.33))
        sigma_efficiency_estimated = np.sqrt(np.diagonal(pcov))
        print("popt: {}".format(popt))

        xplot = np.linspace(20, 1600, 1000)

        plt.plot(xplot,self.efficiency_estimated_popt(xplot,*popt),'r-', color='red')
        #plt.plot(xplot,self.efficiency_estimated_popt(xplot,*(popt+sigma_efficiency_estimated)), color='blue', linewidth=0.4, linestyle='--')
        #plt.plot(xplot,self.efficiency_estimated_popt(xplot,*(popt-sigma_efficiency_estimated)), color='green', linewidth=0.4, linestyle='--')

        #full_width = np.abs(self.efficiency_estimated(xplot,*(popt+sigma_efficiency_estimated))-self.efficiency_estimated(xplot,*(popt-sigma_efficiency_estimated))) #full width of confidence band
        #percent_uncert = full_width/(2*self.efficiency_estimated(xplot,*popt))
        #print("Percentage uncertainty of efficiency {}".format(percent_uncert))

        ##### OBS if want squares and not circles as uncertainty points, use fmt="rs--"
        #plt.errorbar(self.E_Cs[-1], self.eps_Cs[-1], color='green', linewidth=0.001,yerr=self.sigma_eps_Cs[-1], elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        plt.errorbar(self.E_Cs, self.eps_Cs, color='green', linewidth=0.001,yerr=self.sigma_eps_Cs[-1], elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        plt.errorbar(self.E_Ba, self.eps_Ba, color='magenta', linewidth=0.001, yerr=self.sigma_eps_Ba, elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        plt.errorbar(self.E_Eu, self.eps_Eu,color='blue', linewidth=0.001, yerr=self.sigma_eps_Eu, elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')

        #plt.plot(self.E_Cs[-1], self.eps_Cs[-1], '.', color='magenta')
        plt.plot(self.E_Cs, self.eps_Cs, '.', color='magenta')
        plt.plot(self.E_Ba, self.eps_Ba,'.', color='green')
        plt.plot(self.E_Eu, self.eps_Eu, '.', color='blue')
        #plt.yscale('log')
        plt.legend(['Cs137', 'Ba133', 'Eu152', 'fit', 'sigma +', 'sigma -'], loc='best')
        #plt.axis((20,1600,0,0.016))

        #chi = np.sum((efficiency_measured-efficiency_estimated)/error)**2 / (len(efficiency_measured-5)) #5 DEGREES OF FREEDOM
        chi = np.sum(((self.eps - self.efficiency_estimated(self.E))/self.sigma_eps)**2 / len(self.eps) -5)
        print("chi: {}".format(chi))
        plt.xlabel('Energy, keV')
        plt.ylabel('Efficiency')
        plt.legend(['Fit','Cs137', 'Ba133', 'Eu152'], loc='best')

        plt.title('Efficiency calibration of det {}'.format(self.detector))
        plt.savefig("Efficiency_curves/{}".format(self.detector), dpi=300)
        plt.show()
        plt.clf()

    #def chi(self):
    #    np.sum(((self.eps - self.efficiency_estimated(self.E))/self.sigma_eps)**2 / len(self.eps) -5)
        #chi= np.sum(((self.eps - self.efficiency_estimated(self.E))/self.sigma_eps)**2 /len(self.eps - 5))#5 DEGREES OF FREEDOM

        #efficiency_measured-efficiency_estimated)/error)**2 / (len(efficiency_measured-5)) #5 DEGREES OF FREEDOM

#list_detector = ['HPGE1_10','HPGE1_30','HPGE2_10','HPGE2_30', 'IDM1_53','IDM2_32','IDM3_40','IDM4_25', 'room131_5', 'room131_10', 'room131_15', 'room131_18', 'room131_22', 'room131_30', 'room131_40', 'room131_50', 'room131_60'] # ['HPGE', 'HPGE2', 'IDM1', 'IDM2', 'IDM3', 'IDM4', 'room_131_5cm', 'room_131_10cm', 'room_131_15cm', 'room_131_18cm', 'room_131_22cm', 'room_131_30cm', 'room_131_40cm', 'room_131_50cm', 'room_131_60cm']
#for i in list_detector:
#    x=Efficiency_calculations(i)  #Question: can a function name
    #print(x)
    #print(x.__dict__)
#    print(x.plot_func())
    #print(x.plot_data())#.efficiency_measured())
x=Efficiency_calculations('HPGE1_10')  #Question: can a function name
print(x.plot_func())
#chi = x.chi()
#print("Chi: {}".format(chi))
