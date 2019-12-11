import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit

#import fme, os
import os


from foil_info import *

m=60; h=m*60; d=h*24; y=d*356  #converters to seconds

"""
NOTES ON PROGRAM
If driving single decay, activate plot in function, and assign reaction function in end of program.
If driving two step with known parent activity, turn of plot in single decay function. Assign parent reaction and daughter reaction.
If driving two step with unknown parent activity, assign reaction function in end  of program.

All info should be in reaction functions, ie Cu_62Zn, however, guesses are assigned in decay functions. Might need to be changed for some reactions?
"""


####FUNCTIONS

def A0_single_decay(filename_activity_time, lambda_, makePlot=False):
    ID = filename_activity_time[-8:-4]
    Nucleus = filename_activity_time[-12:-8]
    print('foil{}_{}'.format(ID, Nucleus))
    time = np.genfromtxt(filename_activity_time, delimiter=',', usecols=[0]) #hours since e.o.b
    A = np.genfromtxt(filename_activity_time, delimiter=',', usecols=[1])
    sigma_A = np.genfromtxt(filename_activity_time, delimiter=',', usecols=[2])

    index = ~(np.isnan(A) | np.isnan(sigma_A))  #if either eps OR sigma eps is NaN

    t = np.max(time[index])

    xplot = np.linspace(0,t,1000)

    A0_guess=600
    #Single decay mode
    def direct_decay(time, A0_guess):
        A_est=A0_guess*np.exp(-lambda_*time)
        return A_est

    popt, pcov=curve_fit(direct_decay, time[index]*3600, A[index], p0=A0_guess, sigma=sigma_A[index], absolute_sigma=True)
    sigma_activity_estimated = np.sqrt(np.diagonal(pcov))   #Uncertainty in the fitting parameters
    full_width = np.abs(direct_decay(xplot*3600,*(popt+sigma_activity_estimated))-direct_decay(xplot*3600,*(popt-sigma_activity_estimated))) #full width of confidence band
    percent_uncert = full_width/(2*direct_decay(xplot*3600,*popt))  #Distance from fitted line. Half width, sigma_A0/A0 along the line.
    #print("Relative uncertainty of A0 (estimated) {}".format(percent_uncert[0]))
    sigma_A0_estimated = (full_width/2)[0]  #uncertainty in the estimated A0. just taking out the first point.
    A0_estimated = direct_decay(0, popt)
    #print("Activity: {} ({})".format(A0_estimated,sigma_A0_estimated))


    if makePlot == True:
        plt.plot(xplot,direct_decay(xplot*3600,*popt),'r-', color='red')
        plt.plot(xplot,direct_decay(xplot*3600,*(popt+sigma_activity_estimated)), color='blue', linewidth=0.4)
        plt.plot(xplot,direct_decay(xplot*3600,*(popt-sigma_activity_estimated)), color='green', linewidth=0.4)
        plt.plot(time[index],A[index], '.')
        plt.errorbar(time[index], A[index], color='green', linewidth=0.001,yerr=sigma_A[index], elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        plt.title('Activity for foil {} nucleus {}'.format(ID, Nucleus) )
        plt.xlabel('time since eob, hours')
        plt.ylabel('Activity, Bq')
        save_results_to = os.getcwd()+'/activity_curves/'
        #np.savetxt("{}.csv".format(save_results_to +  reaction), np.array((A0, sigma_A0)), delimiter=",")
        plt.savefig('{}foil_{}.png'.format(save_results_to, Nucleus+ID), dpi=300)
        plt.show()




    return A0_estimated, sigma_A0_estimated

def A0_double_decay_unknown_parent(filename_activity_time, lambda_parent, lambda_daughter, makePlot=False):
    ID = filename_activity_time[-8:-4]
    Nucleus = filename_activity_time[-12:-8]
    print('foil{}_{}'.format(ID, Nucleus))

    time = np.genfromtxt(filename_activity_time, delimiter=',', usecols=[0]) #hours since e.o.b
    A = np.genfromtxt(filename_activity_time, delimiter=',', usecols=[1])
    sigma_A = np.genfromtxt(filename_activity_time, delimiter=',', usecols=[2])

    index = ~(np.isnan(A) | np.isnan(sigma_A))  #if either eps OR sigma eps is NaN
    t = np.max(time[index])
    xplot = np.linspace(0,t,1000)

    #non direct decay
    A0_parent_guess=1600; A0_daughter_guess=80000

    def non_direct_decay(time, A0_parent_guess, A0_daughter_guess):  #if there are no gammas to be detected from parent
        A_est = A0_parent_guess*lambda_daughter / (lambda_parent-lambda_daughter) *( np.exp(-lambda_daughter*time)-np.exp(-lambda_parent*time)) + A0_daughter_guess*np.exp(-lambda_daughter *time)
        return A_est

    popt, pcov = curve_fit(non_direct_decay, time[index]*3600, A[index], p0=np.array((A0_parent_guess, A0_daughter_guess)), sigma=sigma_A[index], absolute_sigma=True)
    sigma_activity_estimated = np.sqrt(np.diagonal(pcov))   #Uncertainty in the fitting parameters
    full_width = np.abs(non_direct_decay(xplot*3600,*(popt+sigma_activity_estimated))-non_direct_decay(xplot*3600,*(popt-sigma_activity_estimated))) #full width of confidence band
    percent_uncert = full_width/(2*non_direct_decay(xplot*3600,*popt))  #Distance from fitted line. Half width, sigma_A0/A0 along the line.
    sigma_A0_estimated_isomer = (full_width/2)[0]
    sigma_A0_estimated_ground_state = (full_width/2)[1]  #uncertainty in the estimated A0. just taking out the first point.
    A0_estimated_isomer = popt[0]
    A0_estimated_ground_state = popt[1]
    sigma_A0_estimated = full_width/2 #for 1: isomer, 2: gs of 58Co
    A0_estimated = non_direct_decay(0, popt[0], popt[1])  #58Co  #for 1: isomer, 2: gs of 58Co


    #plot
    if makePlot == True:
        plt.plot(xplot, non_direct_decay(xplot*3600,*popt), 'r-', color='red')
        plt.plot(xplot,non_direct_decay(xplot*3600,*(popt+sigma_activity_estimated)), color='blue', linewidth=0.4)
        plt.plot(xplot,non_direct_decay(xplot*3600,*(popt-sigma_activity_estimated)), color='green', linewidth=0.4)
        plt.plot(time[index],A[index], '.')
        plt.errorbar(time[index], A[index], color='green', linewidth=0.001,yerr=sigma_A[index], elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        plt.title('Activity for foil {} nucleus {}'.format(ID, Nucleus) )
        plt.xlabel('time since eob, hours')
        plt.ylabel('Activity, Bq')
        save_results_to = os.getcwd()+'/activity_curves/'
        #np.savetxt("{}.csv".format(save_results_to +  reaction), np.array((A0, sigma_A0)), delimiter=",")
        plt.savefig('{}foil_{}.png'.format(save_results_to, Nucleus+ID), dpi=300)
        #plt.savefig('foil_{}_{}'.format(Nucleus, ID), dpi=300)
        plt.show()

    #return (A0_estimated)
    return A0_estimated_isomer, sigma_A0_estimated_isomer, A0_estimated_ground_state, sigma_A0_estimated_ground_state

def A0_double_decay_known_parent(filename_activity_time, A0_parent, lambda_parent, lambda_daughter, makePlot=False):
    ID = filename_activity_time[-8:-4]
    Nucleus = filename_activity_time[-12:-8]
    print('foil{}_{}'.format(ID, Nucleus))

    time = np.genfromtxt(filename_activity_time, delimiter=',', usecols=[0]) #hours since e.o.b
    A = np.genfromtxt(filename_activity_time, delimiter=',', usecols=[1])
    sigma_A = np.genfromtxt(filename_activity_time, delimiter=',', usecols=[2])

    index = ~(np.isnan(A) | np.isnan(sigma_A))  #if either eps OR sigma eps is NaN
    t = np.max(time[index])
    xplot = np.linspace(0,t,1000)

    #non direct decay
    A0_daughter_guess=80000

    def non_direct_decay_known_parent(time,A0_daughter_guess):  #if there are no gammas to be detected from parent
        A_est = A0_parent *np.exp(-lambda_parent*time) *lambda_daughter / (lambda_daughter-lambda_parent) * (np.exp(-lambda_daughter*time)-np.exp(-lambda_parent*time)) + A0_daughter_guess*np.exp(-lambda_daughter*time)
        #A_est = A0_parent_guess*lambda_daughter / (lambda_parent-lambda_daughter) *( np.exp(-lambda_daughter*time)-np.exp(-lambda_parent*time)) + A0_daughter_guess*np.exp(-lambda_daughter *time)
        return A_est

    popt, pcov=curve_fit(non_direct_decay_known_parent, time[index]*3600, A[index], p0=A0_daughter_guess, sigma=sigma_A[index], absolute_sigma=True)
    sigma_activity_estimated = np.sqrt(np.diagonal(pcov))   #Uncertainty in the fitting parameters
    full_width = np.abs(non_direct_decay_known_parent(xplot*3600,*(popt+sigma_activity_estimated))-non_direct_decay_known_parent(xplot*3600,*(popt-sigma_activity_estimated))) #full width of confidence band
    percent_uncert = full_width/(2*non_direct_decay_known_parent(xplot*3600,*popt))  #Distance from fitted line. Half width, sigma_A0/A0 along the line.
    #print("Relative uncertainty of A0 (estimated) {}".format(percent_uncert[0]))
    sigma_A0_estimated = (full_width/2)[0]  #uncertainty in the estimated A0. just taking out the first point.
    A0_estimated = non_direct_decay_known_parent(0, popt)
    #print("Activity: {} ({})".format(A0_estimated,sigma_A0_estimated))



    #plot
    if makePlot == True:
        plt.plot(xplot, non_direct_decay_known_parent(xplot*3600,*popt), 'r-', color='red')
        plt.plot(xplot,non_direct_decay_known_parent(xplot*3600,*(popt+sigma_activity_estimated)), color='blue', linewidth=0.4)
        plt.plot(xplot,non_direct_decay_known_parent(xplot*3600,*(popt-sigma_activity_estimated)), color='green', linewidth=0.4)
        plt.plot(time[index],A[index], '.')
        plt.errorbar(time[index], A[index], color='green', linewidth=0.001,yerr=sigma_A[index], elinewidth=0.5, ecolor='k', capthick=0.5)   # cap thickness for error bar color='blue')
        plt.title('Activity for foil {} nucleus {}'.format(ID, Nucleus) )
        plt.xlabel('time since eob, hours')
        plt.ylabel('Activity, Bq')
        save_results_to = os.getcwd()+'/activity_curves/'
        #np.savetxt("{}.csv".format(save_results_to +  reaction), np.array((A0, sigma_A0)), delimiter=",")
        plt.savefig('{}foil_{}.png'.format(save_results_to, Nucleus+ID), dpi=300)
        plt.show()

    #return (A0_estimated)
    return A0_estimated, sigma_A0_estimated




def single_decay_data(func, reaction, n, Save_csv=False):  #function, string
    list, lambda_ = func
    A0 = np.zeros(n); sigma_A0 = np.zeros(n)
    for i,e in enumerate(list):
        A0_estimated, sigma_A0_estimated = A0_single_decay(e, lambda_)
        A0[i] = A0_estimated; sigma_A0[i] = sigma_A0_estimated
        if Save_csv == True:
            save_results_to = os.getcwd()+'/activity_csv/'
            np.savetxt("{}.csv".format(save_results_to +  reaction), np.array((A0, sigma_A0)), delimiter=",")
        print("A0: {}, sigmaA0: {}".format(A0_estimated, sigma_A0_estimated))
        print("*****************************************")
def two_step_kp_data(func_parent, func_daughter, reaction, n, Save_csv=False):
    list_parent, lambda_parent = func_parent
    A0 = np.zeros(n); sigma_A0 = np.zeros(n)
    A0_list = []
    for i,e in enumerate(list_parent):
        A0_estimated_parent, sigma_A0_estimated_parent = A0_single_decay(e,lambda_parent, makePlot=False)
        A0_list.append(A0_estimated_parent)  #add is all A0's for Ni56 from single decay function.
    list_daughter, lambda_parent, lambda_daughter = func_daughter
    for i,e in enumerate(list_daughter):
        A0_estimated_daughter, sigma_A0_estimated_daughter = A0_double_decay_known_parent(e, A0_list[i], lambda_parent, lambda_daughter, makePlot=True)
        A0[i] = A0_estimated_daughter; sigma_A0[i] = sigma_A0_estimated_daughter
        if Save_csv == True:
            save_results_to = os.getcwd()+'/activity_csv/'
            np.savetxt("{}.csv".format(save_results_to +  reaction), np.array((A0, sigma_A0)), delimiter=",")
        #print("A0: {}, sigmaA0: {}".format(A0_estimated_daughter, sigma_A0_estimated_daughter))
        #print("*****************************************")
    print(A0)
def two_step_up_data(func, reaction_parent, reaction_daughter, n, Save_csv=False):
    list, lambda_parent, lambda_daughter = func
    A0_parent = np.zeros(n); sigma_A0_parent = np.zeros(n)
    A0_daughter = np.zeros(n); sigma_A0_daughter = np.zeros(n)
    for i, e in enumerate(list):
        A0_estimated_parent, sigma_A0_estimated_parent, A0_estimated_daughter, sigma_A0_estimated_daughter = A0_double_decay_unknown_parent(e, lambda_parent, lambda_daughter, makePlot=True)
        A0_parent[i] = A0_estimated_parent
        sigma_A0_parent[i] = sigma_A0_estimated_parent
        A0_daughter[i] = A0_estimated_daughter
        sigma_A0_daughter[i] = sigma_A0_estimated_daughter
        if Save_csv == True:
            save_results_to = os.getcwd()+'/activity_csv/'
            np.savetxt("{}.csv".format(save_results_to +  reaction_parent), np.array((A0_parent, sigma_A0_parent)), delimiter=",")
            np.savetxt("{}.csv".format(save_results_to +  reaction_daughter), np.array((A0_daughter, sigma_A0_daughter)), delimiter=",")
            #np.savetxt("{}.csv".format(reaction_parent), np.array((A0_parent, sigma_A0_parent)), delimiter=",")
            #np.savetxt("{}.csv".format(reaction_daughter), np.array((A0_daughter, sigma_A0_daughter)), delimiter=",")
        print("Isomer -  A0: {}, sigma A0: {}".format(A0_estimated_parent, sigma_A0_estimated_parent))
        print("Ground state -  A0: {}, sigma A0: {}".format(A0_estimated_daughter, sigma_A0_estimated_daughter))
        print("***************************************")





#single_decay_data(Cu_62Zn(), "Cu_62Zn", 10, Save_csv=True)
#single_decay_data(Cu_63Zn(), "Cu_63Zn", 10, Save_csv=True)
#single_decay_data(Cu_65Zn(), "Cu_65Zn", 10, Save_csv=True)
#single_decay_data(Ni_57Co(), "Ni_57Co", 10, Save_csv=True)
#single_decay_data(Ni_61Cu(), "Ni_61Cu", 10, Save_csv=True)
#single_decay_data(Cu_65Zn(), "Cu_65Zn", 10, Save_csv=True)
two_step_kp_data(Ni_56Ni(), Ni_56Co(), "Ni_56Co", 10, Save_csv= True)
two_step_up_data(Ni_58Co(), "Ni_58mCo", "Ni_58Co", 10, Save_csv = True)
#two_step_up_data(Cu_52Mn(), "Cu_52mMn", "Ni_52Mn", 10)
#single_decay_data(Fe_56Co(), "Fe_56Co", 3, Save_csv=True)
