from jan20_CrossSections import CrossSections
from ZieglerFiles import ziegler_files
from single_decay_A0 import *
from des19_BeamCurrent import *
#from foil_info import *
#from beam_current_FoilReact import *
#from ZieglerFiles import ziegler_files

files,names = ziegler_files()
CS = CrossSections(files[57])
#path = os.getcwd() + '/activity_csv/'
#file = path + 'Ni_61Cu.csv'

#file = path + 'Cu_57Ni.csv/'
#x = CS.cross_section(Cu_57Ni(), 'Cu', 'Cu_57Ni.csv', 10)

#x = CS.cross_section(Ir_193mPt(), 'Ir', 'Ir_193mPt.csv', 10, 'Ir_193mPt', plot_CS=True)
#x = CS.cross_section(Cu_57Ni(), 'Cu', 'Cu_57Ni.csv', 10, 'Cu_57Ni', plot_CS=True)
#x = CS.cross_section(Cu_64Cu(), 'Cu', 'Cu_64Cu.csv', 10, 'Cu_64Cu', plot_CS=True)


#single_decay_data(Cu_64Cu(), "Cu_64Cu", 10, Save_csv=True)

def get_vals(react_func, target, csv_file, n, reaction):
    CS_class = CrossSections(files[25])
    A0, E, CS, I = CS_class.cross_section(react_func, target, csv_file, n, reaction, plot_CS=True)
    print("A0:", A0)
    print("E:", E)
    print("CS:", CS)
    print("I:", I)

#get_vals(Cu_64Cu(), 'Cu', 'Cu_64Cu.csv', 10, 'Cu_64Cu')
#get_vals(Cu_65Ni(), 'Ni', 'Cu_65Ni.csv', 10, 'Cu_65Ni')
