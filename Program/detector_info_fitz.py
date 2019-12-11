import numpy as np, matplotlib.pyplot as plt
import sys
from dateutil import parser
#sys.path.insert(0, '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Calibration/Calibration_files')
#sys.path.insert(0,'/Users/hannah/Documents/UIO/Masteroppgaven/Data/Calibration/room131/')

#path = 'Master_DataAnalysis/Calibration/fitz_cal_files/'

#f_Cs137=path + 'HPGE1/AO031819_Cs137_10cm_HPGE.txt'
#print(f_Cs137)


#assigning each detector with cal files, live time, real time, and dates of spectra

class Detector_Information_fitz:

    def __init__(self):

        self.date_cal="01/01/2009 12:00:00" #date of source calibrated, same for all sources
        #path to fitz calibration files
        self.path = '/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/'

    def calc_time_delay(self, date_spec):

        from dateutil import parser
        #date_cal = '01-Jan-2009 12:00'
        date_cal="01/01/2009 12:00:00" #date of source calibrated, same for all sources
        dt_cal = parser.parse(date_cal)

        #self.date_spec=[self.date_spec_Cs, self.date_spec_Ba, self.date_spec_Eu]
        dt_spec=[]
        time_d=[]
        time_delay=[]  #time delay after calibration of sources, Cs, Ba, Eu, Na

        for i in date_spec:
        	dt_spec.append(parser.parse(i))
        for j in dt_spec:
        	time_d.append(j-dt_cal)
        for k in time_d:
        	time_delay.append(k.total_seconds())

        self.time_delay=np.array(time_delay)

        return self.time_delay

    def extracting_from_fitz(self, files):
        live_time = []
        date_spec = []
        for i in files:
            f=open(i, 'r')
            lines = f.readlines()
            livetime = lines[21]
            datespec = lines[20]
            lt = livetime.split()
            lt = np.float(lt[-1])
            dt = datespec.split()
            date=dt[1]
            time=dt[-1]
            list=[date,time]
            date_time=','.join(list)
            date_time = date_time.replace(',', ' ')
            date_time = parser.parse(date_time)
            date_spec = []
            live_time.append(lt)
            for j in files:
                datetime=str(date_time)
                date_spec.append(datetime)

            #live_time=np.array(live_time)
        live_time = np.array((live_time))

        f.close()
        return live_time, date_spec



    def HPGE1_10(self):  #Funker ish
        f_Cs137=self.path + 'HPGE1/AO031819_Cs137_10cm_HPGE.txt'
        f_Ba133=self.path + 'HPGE1/AN031819_Ba133_10cm_HPGE.txt'
        f_Eu152 =self.path + 'HPGE1/AM031819_Eu152_10cm_HPGE.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e5,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        B = np.array((25.878e1,25.5e0,0.21637e-1,50.29e-5, 30.33e0))   #b1=2.78e1, b2=2e0
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works

        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152


    def HPGE1_30(self):  #Funker ish
        f_Cs137 = self.path + 'HPGE1/AD022619_Cs137_31.6cm_HPGE.txt'
        f_Ba133 = self.path + 'HPGE1/AE022619_Ba133_31.6cm_HPGE.txt'
        f_Eu152 = self.path + 'HPGE1/AF022619_Eu152_31.6_HPGE.txt'
        #f_Cs137= '/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/HPGE1/AD022619_Cs137_31.6cm_HPGE.txt'
        #f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/HPGE1/AE022619_Ba133_31.6cm_HPGE.txt'
        #f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/HPGE1/AF022619_Eu152_31.6_HPGE.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        B=np.array((7.18e1,2.50e1,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))

        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def HPGE2_10(self): #Funker
        f_Cs137 = self.path + 'HPGE2/AK031719_Cs137_10cm_HPGE2.txt' #'/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/HPGE2/AK031719_Cs137_10cm_HPGE2.txt'
        f_Ba133 = self.path + 'HPGE2/AI031719_Ba133_10cm_HPGE2.txt' #'/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/HPGE2/AI031719_Ba133_10cm_HPGE2.txt'
        f_Eu152 = self.path + 'HPGE2/AJ031719_Eu152_10cm_HPGE2.txt' #'/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/HPGE2/AJ031719_Eu152_10cm_HPGE2.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-1, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))

        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def HPGE2_30(self):  #Something close :)!!!!!!!!!
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/HPGE2/AF022619_Cs137_32cm_hpge2.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/HPGE2/AD022619_Ba133_32cm_hpge2.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/HPGE2/AE022619_Eu152_32cm_hpge2.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        #B = np.array((0.19967883159005945e-1, 0.8172501496981269e0,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def IDM1_53(self):   #Almost work


        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM1/AE20190307_Cs137_53cm_detIDM1.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM1/AF20190329_Ba133_53cm_detIDM1.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM1/AD20190307_Eu152_53cm_detIDM1.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e5, 0.5172501496981269e-2,0.23290419673636917e-3,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        B = np.array((0.078, 2.0, 0.1637, 5.29e-5, 2.33))

        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def IDM2_32(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM2/AC20190307_Cs137_32cm_IDM2.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM2/AE20190308_Ba133_32cm_IDM2.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM2/AD20190307_Eu152_32cm_IDM2.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-6, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works

        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def IDM3_40(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM3/AF20190308_Cs137_38cm_detIDM3.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM3/AD20190308_Ba133_38cm_detIDM3.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM3/AE20190308_Eu152_38cm_detIDM3.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works

        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def IDM4_25(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM4/AE20190307_Cs137_18cm_detIDM4.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM4/AD20190307_Ba133_18cm_detIDM4.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/IDM4/AF20190308_Eu152_18cm_detIDM4.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))

        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152


    def room131_5(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AK20190201_Cs137_05cm_room131.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AL20190201_Ba133_05cm_room131.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AM20190201_Eu152_05cm_room131.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]
        time_delay = self.calc_time_delay(date_spec)
        #B = np.array((2.78e-1,2.0,0.1837e-1,5.29e-5,2.33))  #WORK ###########, chi = 54732.71384855529
        #B = np.array((9.77788899e-5,6.99999999e-1,0.28399999e-1,5.8999999e-5,2.33e-4))  ######not working, chi = 12905.663007807276
        B = np.array((2.7790000e-2,1.99999999e-1,0.28399999e-1,5.3999999e-5,2.33))  ######### working, chi= 24000 something
        #B = np.array((2.7790000e1,1.99999999e0,0.28399999e0,5.3999999e0,2.33e0))

        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))

        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def room131_10(self):  #Not working
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AO20190201_Cs137_10cm_room131.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AP20190201_Ba133_10cm_room131.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AQ20190201_Eu152_10cm_room131.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]
            #live_time=np.array(live_time)
        live_time = np.array((live_time))


        time_delay = self.calc_time_delay(date_spec)
        print(time_delay)
        B = np.array((2.78e3,0.1e-1,0.1637e-1,5.29e-5,2.33))  #NOT WORK
        B = np.array((12e3,0.9,0.1637e-2,5.29e-4,2.33))
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works

        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152


    def room131_15(self):  #Not working
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AS20190201_Cs137_15cm_room131.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AT20190201_Ba133_15cm_room131.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AU20190201_Eu152_15cm_room131.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((1.5e5,10e-7,0.1637e-1,5.29e-5,2.33e-5))


        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-7,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945e3, 0.9172501496981269e-4,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        #B = np.array((7e-4,0.9e-1,0.1637e-2,5.29e-4,2.33))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152


    def room131_18(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AW20190204_Cs137_18cm_room131.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AX20190204_Ba133_18cm_room131.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/AY20190204_Eu152_18cm_room131.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        B = np.array((12e2,0.9,0.1637e-2,5.29e-4,2.33))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))

        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def room131_22(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BA20190204_Cs137_22cm_room131.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BB20190205_Ba133_22cm_room131.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BC20190205_Eu152_22cm_room131.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        B = np.array((12e4,0.9,0.1637e-2,5.29e-4,2.33))   #B1=9.95, B3=5.29e-3
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))

        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def room131_30(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BD20190205_Cs137_30cm_room131.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BE20190206_Ba133_30cm_room131.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BF20190205_Eu152_30cm_room131.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        B = np.array((2.78e2,2.0e-1,0.1637e-1,5.29e-5,2.33))  #FUNKER IKKE
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        B = np.array((50e-3, 2.5e-2,8.11e-2,5.29e-4,20.8e1))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works

        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

    def room131_40(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BG20190206_Cs137_40cm_room131.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BH20190207_Ba133_40cm_room131.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BI20190208_Eu152_40cm_room131.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        B = np.array((2.78e2,2.0e-1,0.1637e-1,5.29e-5,2.33)) #FUNKER IKKE
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        B = np.array((50e-5, 2.5e-3,8.11e-2,5.29e-4,20.8e1))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152


    def room131_50(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BJ20190206_Cs137_50cm_room131.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BK20190208_Ba133_50cm_room131.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BL20190207_Eu152_50cm_room131.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        B = np.array((350,20e-1,0.11e-2,5.29e-5,2.8))
        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))
        B = np.array((50e-1, 2.5e-3,8.11e-2,5.29e-4,20.8e1))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works
        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152


    def room131_60(self):
        f_Cs137='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BM20190206_Cs137_60cm_room131.txt'
        f_Ba133='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BN20190205_Ba133_60cm_room131.txt'
        f_Eu152 ='/Users/hannah/Documents/UIO/Masteroppgaven/Master_DataAnalysis/Calibration/fitz_cal_files/room131/BO20190206_Eu152_60cm_room131.txt'

        files = [f_Cs137, f_Ba133, f_Eu152]

        live_time = self.extracting_from_fitz(files)[0]
        date_spec = self.extracting_from_fitz(files)[1]

        time_delay = self.calc_time_delay(date_spec)
        #B=np.array((7.18e0,2.50e0,0.1637e-2,5.29e-5,2.33e0))
        #B = np.array((0.19967883159005945e-4, 0.9172501496981269e-1,0.23290419673636917e-1,3.4780001711411926e-5,2.4817855558702493 ))
        #B = np.array((1, 2e-3, 5, 5, 0))

        B = np.array((1e4, 2, 0.11e-1,5.29e-2, 2.8e1))
        B = np.array((50e4, 2.5e-3,8.11e-2,5.29e-4,20.8e1))
        B = np.array((0.078, 2.15e0, 0.2637, 5.29e-5, 2.33))    #Works

        #B = np.array((0.19967883159005945, 0.9172501496981269,0.23290419673636917,3.4780001711411926e-05,2.4817855558702493 ))


        return live_time, time_delay, B, f_Cs137, f_Ba133, f_Eu152

x=Detector_Information_fitz()
#det = x.room131_18()
#print(x.calc_time_delay("01/30/2019 15:58:37"))
#print(x.HPGE1_32)
#print(det)
