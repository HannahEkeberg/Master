import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize_scalar
from scipy import asarray as ar,exp
from scipy.stats import norm

path = '/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Beamspot_film/'
filename = path +'ss2_beamprofile.csv'
textimage = path + 'ss1_textimage.txt'



E = np.genfromtxt(filename, delimiter=',', usecols=0, skip_header=1)
I = np.genfromtxt(filename, delimiter=',', usecols=1, skip_header=1)

def gauss_fit(E,I,val_min, val_max):
    #n = int((E_peak+peak_range)-(E_peak-peak_range))
    counts = []#np.zeros(int(n))
    E = []#np.zeros(int(n))

    #val_max = E_peak+peak_range
    #val_min = E_peak-peak_range
    index = []
    for i, v in enumerate(energy): #index, value
        if v >= val_min and v <= val_max:
            index.append(i)
            E.append(energy[i])
            counts.append(spectrum[i])
    counts=np.array((counts)); E=np.array((E))

    n = len(E)
    #mean = sum(E, counts)/n
    #sigma = sum(counts*(E-mean)**2)/n
    #print(mean, sigma)
    (mu,sigma)=norm.fit(E)
    #print(mu, sigma)
    #X = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)

    def gauss(x,a, x0, sigma):
        #return 1/(np.sqrt(2*np.pi)) * 1/sigma * np.exp(-(x-mu)**2 /(2*sigma**2))
        return a*np.exp(-(x-x0)**2 / (2*sigma**2))

    popt, pcov = curve_fit(gauss, E, counts, p0=[max(counts), mu, sigma], absolute_sigma=True)
    fm = lambda E: -gauss(E, *popt)
    r = minimize_scalar(fm, bounds=(min(E), max(E)))
    #print(r["x"], gauss(r["x"], *popt))
