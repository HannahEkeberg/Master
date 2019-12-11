import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit

m=60; h=m*60; d=h*24; y=d*356  #converters to seconds

fnames = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/fnames.txt'
with open (fnames, "r") as myfile:
    data=myfile.readlines()
    #print(type(data))


path = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/'


####Cupper foils########
def Cu_62Zn(): #check,  #mon, single
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_129.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_229.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_329.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_429.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_529.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_629.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_729.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_829.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_929.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_1029.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(9.193*h)
    return list, lambda_    #mon, single

def Cu_63Zn(): #check #mon, single
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_129.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_229.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_329.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_429.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_529.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_629.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_729.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_829.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_929.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/63Zn_1029.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(38.47*m)
    return list, lambda_    #mon, single

def Cu_65Zn(): #check,  #mon, single
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_129.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_229.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_329.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_429.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_529.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_629.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_729.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_829.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_929.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/65Zn_1029.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(243.93*d)
    return list, lambda_    #mon

def Cu_52Mn(): #non-mon, BUT WAS PRODUCED, needs work ?????????
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_129.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_229.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_329.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_429.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_529.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_629.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_729.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_829.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_929.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/52Mn_1029.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_parent = np.log(2)/(21.1*m) #52mMn
    lambda_daughter = np.log(2)/(5.591*d) #52Mn
    return list, lambda_parent, lambda_daughter


def Cu_56Co():
    foil1 = path + '56Co_129.dat'
    foil2 = path + '56Co_229.dat'
    foil2 = path + '56Co_329.dat'
    foil2 = path + '56Co_429.dat'
    foil2 = path + '56Co_529.dat'
    foil2 = path + '56Co_629.dat'
    foil2 = path + '56Co_729.dat'
    foil2 = path + '56Co_829.dat'
    foil2 = path + '56Co_929.dat'
    foil2 = path + '56Co_1029.dat'



    #pass

def Cu_57Co():
    pass
def Cu_57Ni():
    pass
def Cu_58Co():
    pass
def Cu_59Fe():
    pass
def Cu_60Co():
    pass
def Cu_61Co():
    pass
def Cu_61Cu():   #in theory, decay from 61Zn but half life in order seconds...
    foil1 = path + '/61Cu_129.dat'
    foil2 = path + '/61Cu_229.dat'
    foil3 = path + '/61Cu_329.dat'
    foil4 = path + '/61Cu_429.dat'
    foil5 = path + '/61Cu_529.dat'
    foil6 = path + '/61Cu_629.dat'
    foil7 = path + '/61Cu_729.dat'
    foil8 = path + '/61Cu_829.dat'
    foil9 = path + '/61Cu_929.dat'
    foil10 = path + '/61Cu_1029.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(3.339*h)
    return list, lambda_

def Cu_64Cu():
    pass
def Cu_65Ni():
    pass

#####NICKEL FOILS######

def Ni_56Ni(): #non-mon, single
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_128.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_228.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_328.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_428.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_528.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_628.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_728.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_828.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_928.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Ni_1028.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(6.075*d)
    return list, lambda_

def Ni_56Co(): #mon, two step, known parent activity 56Ni
    #type = "tskp"
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_128.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_228.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_328.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_428.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_528.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_628.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_728.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_828.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_928.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_1028.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_daughter = np.log(2)/(77.236*d) #56Co
    lambda_parent = np.log(2)/(6.075*d) #56Ni
    return list, lambda_parent, lambda_daughter

def Ni_58Co(): #mon, two step, unkown parent activity 58mCo
    #type = "tsup"
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_128.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_228.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_328.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_428.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_528.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_628.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_728.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_828.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_928.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/58Co_1028.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_isomer = np.log(2)/(9.10*h)
    lambda_ground_state = np.log(2)/(77.236*d)
    return list, lambda_isomer, lambda_ground_state

def Ni_61Cu(): #mon, single
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_128.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_228.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_328.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_428.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_528.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_628.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_728.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_828.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_928.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/61Cu_1028.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(3.339*h)
    return list, lambda_

def Ni_57Co(): #single
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_128.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_228.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_328.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_428.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_528.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_628.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_728.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_828.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_928.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/57Co_1028.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(271.74*d)
    return list, lambda_

def Ni_57Ni():
    pass

def Ni_55Co(): #single, beta- from 55Ni too short half life..
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_128.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_228.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_328.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_428.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_528.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_628.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_728.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_828.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_928.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/55Co_1028.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(17.53*h)
    return list, lambda_   #gives weird activities

def Ni_52mMn():
    pass

def Ni_52Mn():
    pass

def Ni_54Mn():
    pass

def Ni_59Fe():
    pass

def Ni_60Cu():
    pass

def Ni_60mCo():
    pass

def Ni_64Cu():
    pass

def Ni_65Ni():
    pass
#####IRON FOILS########

def Fe_56Co(): #mon, Single
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_126.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_226.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/56Co_326.dat'
    list = [foil1, foil2, foil3]
    lambda_ = np.log(2)/(77.236*d)
    return list, lambda_

def Fe_48V():
    pass

def Fe_51Cr():
    pass
def Fe_51Mn():
    pass
def Fe_52mMn():
    pass
def Fe_52Mn():
    pass
def Fe_53Fe():
    pass
def Fe_54Mn():
    pass
def Fe_55Co():
    pass
def Fe_57Co():
    pass
def Fe_58Co():
    pass
def Fe_59Fe():
    pass

###Iridium foils####

def Ir_183Ta():
    pass
def Ir_186Re():
    pass
def Ir_186Ta():
    pass
def Ir_187W():
    pass
def Ir_188Ir():
    pass
def Ir_188mRe():
    pass
def Ir_188Pt():
    pass
def Ir_188Re():
    pass
def Ir_189Ir():
    pass
def Ir_189Pt():
    pass
def Ir_189Re():
    pass
def Ir_189W():
    pass
def Ir_190Ir():
    pass
def Ir_190mRe():
    pass
def Ir_190Re():
    pass
def Ir_191Pt():
    pass
def Ir_192Ir():
    pass
def Ir_193mPt(): #HELLUUUUU
    foil1 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_177.dat'
    foil2 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_277.dat'
    foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_377.dat'
    foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_477.dat'
    foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_577.dat'
    foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_677.dat'
    foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_777.dat'
    foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_877.dat'
    foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_977.dat'
    foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/193mPt_1077.dat'
    list = [foil1, foil2, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(4.33*d)  #days to seconds
    return list, lambda_    #mon, single

def Ir_194Ir():
    pass
def Ir_194m2Ir():
    pass
