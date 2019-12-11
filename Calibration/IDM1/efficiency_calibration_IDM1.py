import matplotlib.pyplot as plt
from numpy import exp, loadtxt, pi, sqrt
import numpy as np
from scipy.optimize import curve_fit


#from lmfit import Model

#First part must be changed depending on which detector

f_Cs137='AI20190131_Cs137_10cm_IDM1.sin'
f_Ba133='AJ20190131_Ba133_10cm_IDM1.sin'
f_Eu152='AL20190204_Eu152_10cm_IDM1.sin'
#f_Na22='AD20190131_Na22_10cm_DetHPGE.sin'

date_spec= ["01/31/2019 16:13:44","01/31/2019 16:23:33","02/04/2019 10:39:00"]
real_time=np.array((384,6217,14103)) #seconds  
live_time=np.array((378,6097,13767)) #seconds 

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#Same for all detectors

t_half=np.array((30.08,10.551,13.517))*365*24*3600 # seconds, Cs137, Ba133, Eu152
t_half_uncert = np.array((0.09,0.011, 0.014))*365*24*3600
A_0=np.array((38.55,39.89,39.29))*1e3 # Source activity Bq, Cs,Ba,Eu
Intensity_perc = np.array((0.8510,0.6205,0.2087)) #Intensity of strongest lines in spectras

from dateutil import parser
date_cal="01/01/2009 12:00:00" #date of source calibrated, same for all sources
dt_cal = parser.parse(date_cal)

dt_spec=[]
time_d=[]
time_delay=[]  #time delay after calibration of sources, Cs, Ba, Eu, Na

for i in date_spec:
	dt_spec.append(parser.parse(i))
for j in dt_spec:
	time_d.append(j-dt_cal)
for k in time_d: 
	time_delay.append(k.total_seconds())

time_delay=np.array(time_delay)



#f_list=[]
#f_list.extend([f_Cs137, f_Ba133, f_Eu152, f_Na22])






def cal(filename, t_half, A_0, live_time, real_time, time_delay, Intensity_perc, t_half_uncert):
	channel=np.loadtxt(filename, usecols=[1], skiprows=1)
	sigma_channel=np.loadtxt(filename, usecols=[2], skiprows=1)
	N_c=np.loadtxt(filename, usecols=[3], skiprows=1) #net area
	sigma_Nc=np.loadtxt(filename, usecols=[4], skiprows=1)
	E=np.loadtxt(filename, usecols=[5], skiprows=1)
	sigma_E=np.loadtxt(filename, usecols=[6], skiprows=1)
	I=np.loadtxt(filename, usecols=[7], skiprows=1) #intensities. 
	sigma_I=np.loadtxt(filename, usecols=[8], skiprows=1)

	lamb=np.log(2)/t_half
	sigma_lamb= lamb*(t_half_uncert/t_half) #np.log(2)/t_half_uncert

	sigma_A_0=A_0*0.03 #general uncertainty of calibration sources


	dead_time=real_time-live_time

	I_gamma=I*Intensity_perc/1e4 #Remember in cal files for each source, highest intenstity was set to 10.000 counts
	sigma_I_gamma = sigma_I*Intensity_perc/1e4
	
	eps=N_c*lamb/(A_0*I_gamma*(1-np.exp(-lamb*live_time))*np.exp(-lamb*time_delay))

	sigma_eps=eps*np.sqrt((sigma_lamb/lamb)**2 + (sigma_Nc/N_c)**2 + (sigma_A_0/A_0)**2 + (sigma_I_gamma/I_gamma)**2  )
	

	#print N_c
	return E,eps, sigma_E, sigma_eps


E_Cs,eps_Cs, sigma_E_Cs, sigma_eps_Cs = cal(f_Cs137, t_half[0], A_0[0], live_time[0], real_time[0], time_delay[0], Intensity_perc[0], t_half_uncert[0])
E_Ba,eps_Ba,sigma_E_Ba, sigma_eps_Ba = cal(f_Ba133, t_half[1], A_0[1], live_time[1], real_time[1], time_delay[1], Intensity_perc[1], t_half_uncert[1])
E_Eu,eps_Eu, sigma_E_Eu, sigma_eps_Eu = cal(f_Eu152, t_half[2], A_0[2], live_time[2], real_time[2], time_delay[2], Intensity_perc[2], t_half_uncert[2])
#E_Na,eps_Na,sigma_E_Na, sigma_eps_Na = cal(f_Na22, t_half[3], A_0[3], live_time[3], real_time[3], time_delay[3], Intensity_perc[3])

sigma_eps_Cs, sigma_eps_Ba, sigma_eps_Eu

plt.plot(E_Cs, eps_Cs,'o', color='magenta')
plt.plot(E_Ba, eps_Ba,'o', color='green')
plt.plot(E_Eu, eps_Eu, 'o', color='blue')
#plt.plot(E_Na, eps_Na, 'o')
#plt.plot(E_Cs,eps_Cs,'-.', linewidth=0.4)
#plt.plot(E_Ba, eps_Ba,'-.',linewidth=0.4)
#plt.plot(E_Eu, eps_Eu, '-.',linewidth=0.4)
#plt.plot(E_Na, eps_Na, '-.',linewidth=0.4)
plt.xlabel('Energy, keV')
plt.ylabel('Efficiency')
plt.legend(['Cs137', 'Ba133', 'Eu152'], loc='best')
plt.title('Calibration of IDM1, Berkeley CA')
#plt.show()


x = np.concatenate((E_Cs, E_Ba, E_Eu))#, np.array((E_Na))), axis=1)
y = np.concatenate((eps_Cs, eps_Ba, eps_Eu))
x_err = np.concatenate((sigma_E_Cs, sigma_E_Ba, sigma_E_Eu))
y_err = np.concatenate((sigma_eps_Cs, sigma_eps_Ba, sigma_eps_Eu))


def efficiency(E_gamma, B_0, B_1, B_2, B_3, B_4):
    return B_0*np.exp(-B_1*E_gamma**B_2)*(1-np.exp(-B_3*E_gamma**B_4))


#print sigma_eps_Cs/eps_Cs
#print sigma_eps_Ba/eps_Ba
#print sigma_eps_Eu/eps_Eu

E_gamma=x
popt, pcov=curve_fit(efficiency, x, y, p0=np.array([2.78,2.0,0.1637,5.29e-5,2.33]),sigma=y_err/y) #John's numbers

#popt, pcov=curve_fit(efficiency, x, y, p0=np.array([1e-4,-10,-2,1e3,-1]),sigma=y_err/y)  #different working versions
#popt, pcov=curve_fit(efficiency, x, y, p0=np.array([1e-4,-10,-3,1e3,-1]),)  WORKS!!!!!!!!!!!!!!!!!!!!!!!!

xplot = np.linspace(20, 1600, 1000)

def chi(efficiency_measured, efficiency_approx, error):  #Calculate chi 
	return np.sum(((efficiency_measured-efficiency_approx)/error)**2) / (len(efficiency_measured)- 5 ) #(divided by degrees of freedom)

print "chi=%.4f" %chi(y, efficiency(x,popt[0], popt[1], popt[2], popt[3], popt[4]), y_err)  # Chi with John's numbers closest to 1, works best


plt.plot(xplot,efficiency(xplot,*popt),'r-', color='red')

##### OBS if want squares and not circles as uncertainty points, use fmt="rs--"
plt.errorbar(E_Cs, eps_Cs, color='green', linewidth=0.001, xerr=sigma_E_Cs, yerr=sigma_eps_Cs, elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
plt.errorbar(E_Ba, eps_Ba, color='magenta', linewidth=0.001, xerr=sigma_E_Ba, yerr=sigma_eps_Ba, elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
plt.errorbar(E_Eu, eps_Eu,color='blue', linewidth=0.001, xerr=sigma_E_Eu, yerr=sigma_eps_Eu, elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')

#plt.yscale('log')
plt.legend(['Cs137', 'Ba133', 'Eu152', 'fit'], loc='best')
plt.savefig('Calibration_IDM1', dpi=300)
#plt.show()


#Testing for Cs137  #####################################################################################
def activity_test(lamb, N_c, eps, I_gamma, t_live):
	return (lamb*N_c) / (eps*I_gamma*(1-exp(-lamb*t_live)))

eff_calibrated= efficiency(661,popt[0], popt[1], popt[2], popt[3], popt[4])   #find efficiency at known energy peak, here Cs 
#print eff_calibrated

lamb=np.log(2)/t_half[0]

#print activity_test(np.log(2)/t_half[0],134392.0,eps_Cs, 0.851, live_time[0] ) #activity with measured efficiency
A0= activity_test(np.log(2)/t_half[0],24044.0,eff_calibrated, 0.851, live_time[0] )  #with efficiency found from calibration curve

print "A0=%.4f "%((A0/np.exp(-lamb*time_delay[0]) ) *1e-3)



